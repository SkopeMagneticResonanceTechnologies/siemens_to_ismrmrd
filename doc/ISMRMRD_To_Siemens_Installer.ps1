####################################################################################
#    Script for installing the required dependencies for 
#    the ISMRM Raw Data Format and Siemens converter on Windows 10.
#    
#    Prerequisites:
#        - Windows 10 (64-bit)
#        - Visual Studio 2017 Community (C/C++) installed
#        - 7zip
#   
#   By default the execution of scripts is not permitted, therefore run the 
#   following command in a PowerShell with administrator rights:
#        Set-ExecutionPolicy RemoteSigned
####################################################################################

###################################################################
####################      Functions   #############################
###################################################################

function download_file($url,$destination) {    
    #Let's set up a webclient for all the files we have to download   
    [Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12    
    $client = New-Object System.Net.WebClient    
    $client.DownloadFile($url,$destination)
}

function unzip($zipPath, $destination){    
    $shell = new-object -com shell.application;    
    $zip = $shell.NameSpace($zipPath);    

    if ((Test-Path -path $destination) -ne $True){        
        New-Item $destination -type directory    
    }    

    foreach($item in $zip.items()){  
      $shell.Namespace($destination).copyhere($item)    }

}

function add_path($pathname) {    
    if ($env:path  -match [regex]::escape($pathname)){
        Write-Host "Path $path already exists"    } 
    else {        
        setx PATH "$env:path;$pathname" -m    
    }
}

function Get-Exports {
<#
FuzzySecurity: https://github.com/FuzzySecurity/PowerShell-Suite

.SYNOPSIS
Get-Exports, fetches DLL exports and optionally provides
C++ wrapper output (idential to ExportsToC++ but without
needing VS and a compiled binary). To do this it reads DLL
bytes into memory and then parses them (no LoadLibraryEx).
Because of this you can parse x32/x64 DLL's regardless of
the bitness of PowerShell.
.DESCRIPTION
Author: Ruben Boonen (@FuzzySec)
License: BSD 3-Clause
Required Dependencies: None
Optional Dependencies: None
.PARAMETER DllPath
Absolute path to DLL.
.PARAMETER CustomDll
Absolute path to output file.
.EXAMPLE
C:\PS> Get-Exports -DllPath C:\Some\Path\here.dll
.EXAMPLE
C:\PS> Get-Exports -DllPath C:\Some\Path\here.dll -ExportsToCpp C:\Some\Out\File.txt
#>
    param (
        [Parameter(Mandatory = $True)]
        [string]$DllPath,
        [Parameter(Mandatory = $False)]
        [string]$ExportsToCpp
    )

    Add-Type -TypeDefinition @" 
    using System;
    using System.Diagnostics;
    using System.Runtime.InteropServices;
    using System.Security.Principal;

    [StructLayout(LayoutKind.Sequential)]
    public struct IMAGE_EXPORT_DIRECTORY
    {
        public UInt32 Characteristics;
        public UInt32 TimeDateStamp;
        public UInt16 MajorVersion;
        public UInt16 MinorVersion;
        public UInt32 Name;
        public UInt32 Base;
        public UInt32 NumberOfFunctions;
        public UInt32 NumberOfNames;
        public UInt32 AddressOfFunctions;
        public UInt32 AddressOfNames;
        public UInt32 AddressOfNameOrdinals;
    }

    [StructLayout(LayoutKind.Sequential)]
    public struct IMAGE_SECTION_HEADER
    {
        public String Name;
        public UInt32 VirtualSize;
        public UInt32 VirtualAddress;
        public UInt32 SizeOfRawData;
        public UInt32 PtrToRawData;
        public UInt32 PtrToRelocations;
        public UInt32 PtrToLineNumbers;
        public UInt16 NumOfRelocations;
        public UInt16 NumOfLines;
        public UInt32 Characteristics;
    }

    public static class Kernel32
    {
        [DllImport("kernel32.dll")]
        public static extern IntPtr LoadLibraryEx(
            String lpFileName,
            IntPtr hReservedNull,
            UInt32 dwFlags);
    }
"@

    # Load the DLL into memory so we can reference it like LoadLibrary
    $FileBytes = [System.IO.File]::ReadAllBytes($DllPath)
    [IntPtr]$HModule = [System.Runtime.InteropServices.Marshal]::AllocHGlobal($FileBytes.Length)
    [System.Runtime.InteropServices.Marshal]::Copy($FileBytes, 0, $HModule, $FileBytes.Length)

    # Some Offsets..
    $PE_Header = [Runtime.InteropServices.Marshal]::ReadInt32($HModule.ToInt64() + 0x3C)
    $Section_Count = [Runtime.InteropServices.Marshal]::ReadInt16($HModule.ToInt64() + $PE_Header + 0x6)
    $Optional_Header_Size = [Runtime.InteropServices.Marshal]::ReadInt16($HModule.ToInt64() + $PE_Header + 0x14)
    $Optional_Header = $HModule.ToInt64() + $PE_Header + 0x18

    # We need some values from the Section table to calculate RVA's
    $Section_Table = $Optional_Header + $Optional_Header_Size
    $SectionArray = @()
    for ($i; $i -lt $Section_Count; $i++) {
        $HashTable = @{
            VirtualSize = [Runtime.InteropServices.Marshal]::ReadInt32($Section_Table + 0x8)
            VirtualAddress = [Runtime.InteropServices.Marshal]::ReadInt32($Section_Table + 0xC)
            PtrToRawData = [Runtime.InteropServices.Marshal]::ReadInt32($Section_Table + 0x14)
        }
        $Object = New-Object PSObject -Property $HashTable
        $SectionArray += $Object

        # Increment $Section_Table offset by Section size
        $Section_Table = $Section_Table + 0x28
    }

    # Helper function for dealing with on-disk PE offsets.
    # Adapted from @mattifestation:
    # https://github.com/mattifestation/PowerShellArsenal/blob/master/Parsers/Get-PE.ps1#L218
    function Convert-RVAToFileOffset($Rva, $SectionHeaders) {
        foreach ($Section in $SectionHeaders) {
            if (($Rva -ge $Section.VirtualAddress) -and
                ($Rva-lt ($Section.VirtualAddress + $Section.VirtualSize))) {
                return [IntPtr] ($Rva - ($Section.VirtualAddress - $Section.PtrToRawData))
            }
        }
        # Pointer did not fall in the address ranges of the section headers
        echo "Mmm, pointer did not fall in the PE range.." 
    }

    # Read Magic UShort to determine x32/x64
    if ([Runtime.InteropServices.Marshal]::ReadInt16($Optional_Header) -eq 0x010B) {
        #echo "`n[?] 32-bit Image!" 
        # IMAGE_DATA_DIRECTORY[0] -> Export
        $Export = $Optional_Header + 0x60
    } else {
        #echo "`n[?] 64-bit Image!" 
        # IMAGE_DATA_DIRECTORY[0] -> Export
        $Export = $Optional_Header + 0x70
    }

    # Convert IMAGE_EXPORT_DIRECTORY[0].VirtualAddress to file offset!
    $ExportRVA = Convert-RVAToFileOffset $([Runtime.InteropServices.Marshal]::ReadInt32($Export)) $SectionArray

    # Cast offset as IMAGE_EXPORT_DIRECTORY
    $OffsetPtr = New-Object System.Intptr -ArgumentList $($HModule.ToInt64() + $ExportRVA)
    $IMAGE_EXPORT_DIRECTORY = New-Object IMAGE_EXPORT_DIRECTORY
    $IMAGE_EXPORT_DIRECTORY = $IMAGE_EXPORT_DIRECTORY.GetType()
    $EXPORT_DIRECTORY_FLAGS = [system.runtime.interopservices.marshal]::PtrToStructure($OffsetPtr, [type]$IMAGE_EXPORT_DIRECTORY)
    # Get equivalent file offsets!
    $ExportFunctionsRVA = Convert-RVAToFileOffset $EXPORT_DIRECTORY_FLAGS.AddressOfFunctions $SectionArray
    $ExportNamesRVA = Convert-RVAToFileOffset $EXPORT_DIRECTORY_FLAGS.AddressOfNames $SectionArray
    $ExportOrdinalsRVA = Convert-RVAToFileOffset $EXPORT_DIRECTORY_FLAGS.AddressOfNameOrdinals $SectionArray

    # Loop exports
    $ExportArray = @()
    for ($i=0; $i -lt $EXPORT_DIRECTORY_FLAGS.NumberOfNames; $i++){
        # Calculate function name RVA
        $FunctionNameRVA = Convert-RVAToFileOffset $([Runtime.InteropServices.Marshal]::ReadInt32($HModule.ToInt64() + $ExportNamesRVA + ($i*4))) $SectionArray
        $HashTable = @{
            FunctionName = [System.Runtime.InteropServices.Marshal]::PtrToStringAnsi($HModule.ToInt64() + $FunctionNameRVA)
            ImageRVA = "0x$("{0:X8}" -f $([Runtime.InteropServices.Marshal]::ReadInt32($HModule.ToInt64() + $ExportFunctionsRVA + ($i*4))))" 
            Ordinal = [Runtime.InteropServices.Marshal]::ReadInt16($HModule.ToInt64() + $ExportOrdinalsRVA + ($i*2)) + $EXPORT_DIRECTORY_FLAGS.Base
        }
        $Object = New-Object PSObject -Property $HashTable
        $ExportArray += $Object
    }

    # Optionally write ExportToC++ wrapper output
    if ($ExportsToCpp) {
        Add-Content $ExportsToCpp "EXPORTS" 
        foreach ($Entry in $ExportArray) {
            Add-Content $ExportsToCpp "$($Entry.FunctionName)" 
        }
    }

    # Free buffer
    [Runtime.InteropServices.Marshal]::FreeHGlobal($HModule)
}

# Invokes a Cmd.exe shell script and updates the environment.
function Invoke-CmdScript 
{
<#
.SYNOPSIS
Executes the specified command script and imports the environment into current
PowerShell instance.

.DESCRIPTION
Running development tools at the command line in PowerShell can be a hassle since
they rely on environment varibles and those are set through batch files. This
cmdlet runs those batch files and imports any set environment variables into
the running PowerShell instance.

.PARAMETER script
The required batch file to run.

.PARAMETER parameters
The optional parameters to pass to the batch file.

.NOTES
The original script is by Lee Holmes. I updated the script to make removing environment variables
work.

.LINK
http://www.leeholmes.com/blog/2006/05/11/nothing-solves-everything-%e2%80%93-powershell-and-other-technologies/
https://github.com/Wintellect/WintellectPowerShell
#>
    param
    (
        [Parameter(Mandatory=$true,
                   Position=0,
                   HelpMessage="Please specify the command script to execute.")]
        [string] $script, 
        [Parameter(Position=1)]
        [string] $parameters="" 
    )  

    # Save off the current environment variables in case there's an issue
    $oldVars = $(Get-ChildItem -Path env:\)
    $tempFile = [IO.Path]::GetTempFileName()  

    try
    {
        ## Store the output of cmd.exe. We also ask cmd.exe to output
        ## the environment table after the batch file completes
        cmd /c " `"$script`" $parameters && set > `"$tempFile`" " 

        if ($LASTEXITCODE -ne 0)
        {
            throw "Error executing CMD.EXE: $LASTEXITCODE" 
        }

        # Before we delete the environment variables get the output into a string
        # array.
        $vars = Get-Content -Path $tempFile

        # Clear out all current environment variables in PowerShell.
        Get-ChildItem -Path env:\ | Foreach-Object { 
                        set-item -force -path "ENV:\$($_.Name)" -value "" 
                    }

        ## Go through the environment variables in the temp file.
        ## For each of them, set the variable in our local environment.
        $vars | Foreach-Object {   
                            if($_ -match "^(.*?)=(.*)$")  
                            { 
                                Set-Content -Path "env:\$($matches[1])" -Value $matches[2]
                            } 
                        }
    }
    catch
    {
        "ERROR: $_" 

        # Any problems, restore the old environment variables.
        $oldVars | ForEach-Object { Set-Item -Force -Path "ENV:\$($_.Name)" -value $_.Value }
    }
    finally
    {
        Remove-Item -Path $tempFile -Force -ErrorAction SilentlyContinue
    }
}
###################################################################
####################      Functions End  ##########################
###################################################################


#Clean up
&cls

Write-Host "**********************************************************************" 
Write-Host "**********************************************************************" 
Write-Host "** ISMRMRD Raw Data Format & SIEMENS Data Converter Installation      " 
Write-Host "**********************************************************************" 
Write-Host "**********************************************************************" 
Write-Host "" 
Write-Host "" 

#Storage location for ISMRMRD and converter source files
$library_location = "C:\MRILibraries" 
$download_location = "C:\MRILibraries\downloads" 
$currentFolder = pwd



###################################################################################
### Preparations
###################################################################################
#region Checks 
#The user must be an administrator
If (-NOT ([Security.Principal.WindowsPrincipal] [Security.Principal.WindowsIdentity]::GetCurrent()).IsInRole([Security.Principal.WindowsBuiltInRole] "Administrator")){
    Write-Warning "You do not have Administrator rights to run this script!`nPlease re-run this script as an Administrator!" 
    Break
}

#Check for 7zip
if (-not (test-path "$env:ProgramFiles\7-Zip\7z.exe")) {
    Write-Warning "$env:ProgramFiles\7-Zip\7z.exe needed." 
    Break
}
set-alias sz "$env:ProgramFiles\7-Zip\7z.exe"  

if ((Test-Path -path "C:\Program Files (x86)\Microsoft Visual Studio\2017\Community\VC\Auxiliary\Build\vcvarsall.bat") -ne $True){    
    Write-Warning "Microsoft Visual Studio 2017 Community needed.`nPlease install Visual Studio first!" 
    Break
}
#endregion

#Let's first check if we have the library folder and if not create it
if ((Test-Path -path $library_location) -ne $True){    
    Write-Host "Library location: " $library_location " not found, creating"    
    New-Item $library_location -type directory
}
else{    
    Write-Host "Library location: " $library_location " found." 
}
Write-Host "" 
Write-Host "" 

#Now check if we have the downlaod folder and if not create it
if ((Test-Path -path $download_location) -ne $True){
    Write-Host "Download location: " $download_location " not found, creating"    
    New-Item $download_location -type directory
}
else{      
    Write-Host "Download location: " $download_location " found." 
    Write-Host "Cleaning up..." 
    Remove-Item $download_location\* -Force -Recurse
}
Write-Host "" 
Write-Host "" 



###################################################################################
### Set environment variables for Visual Studio
###################################################################################
Invoke-CmdScript "C:\Program Files (x86)\Microsoft Visual Studio\2017\Community\VC\Auxiliary\Build\vcvarsall.bat" "x64" 
Write-Host "" 
Write-Host "" 



###################################################################################
### Install dependencies
###################################################################################
#Download and install CMAKE
if ((Test-Path -path "C:\Program Files (x86)\CMake") -ne $True){ 
    Write-Host "Downloading cmake (v3.11.1)... (Remember to add cmake to the system PATH of the current user)" 
    download_file "https://cmake.org/files/v3.11/cmake-3.11.1-win32-x86.msi" ($download_location + "\cmake-3.11.1-win32-x86.msi")
    Write-Host "Installing cmake" 
    Start-Process msiexec.exe -Wait -ArgumentList ("/I " + $download_location + "\cmake-3.11.1-win32-x86.msi")
}
else{
    Write-Host "cmake seems to be installed already." 
}
Write-Host "" 
Write-Host "" 

#Download and install Git
if ((Test-Path -path "C:\Program Files\Git") -ne $True){ 
    Write-Host "Downloading Git for Windows (2.17.0-64-bit)... (Use the proposed default settings)" 
    download_file "https://github.com/git-for-windows/git/releases/download/v2.17.0.windows.1/Git-2.17.0-64-bit.exe" ($download_location + "\Git-2.17.0-64-bit.exe")
    Write-Host "Installing Git" 
    Start-Process ($download_location + "\Git-2.17.0-64-bit.exe") -Wait 
}
else{
    Write-Host "Git seems to be installed already." 
}
Write-Host "" 
Write-Host "" 

#Download, unzip, and install HDF5
if ((Test-Path -path "C:\Program Files\HDF Group\HDF5\1.8.9") -ne $True){ 
Write-Host "Downloading HDF5 (1.8.9)... (Remember to add HDF5 to the system PATH of the current user. Do not worry about 'Path too long' warning if you used the default install path.)" 
download_file "https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.9/bin/windows/HDF5189-win64-vs10-shared.zip" ($download_location + "\HDF5189-win64-vs10-shared.zip")
Write-Host "Unzipping HDF5..." 
unzip ($download_location + "\HDF5189-win64-vs10-shared.zip")  "$download_location\hdf5_binaries" 
Write-Host "Installing HDF5" 
Start-Process ($download_location + "\hdf5_binaries\HDF5-1.8.9-win64.exe") -Wait 
}
else{
    Write-Host "HDF5 seems to be installed already." 
}
Write-Host "" 
Write-Host "" 

#FFTW
if ((Test-Path -path "$library_location\fftw3") -ne $True){ 
    Write-Host "Downloading FFTW (3.3.5-dll64)..." 
    download_file "ftp://ftp.fftw.org/pub/fftw/fftw-3.3.5-dll64.zip" ($download_location + "\fftw-3.3.5-dll64.zip")
    Write-Host "Unzipping BOOST..." 
    unzip ($download_location + "\fftw-3.3.5-dll64.zip")  "$library_location\fftw3" 
    cd "$library_location\fftw3" 
    & lib "/machine:X64" "/def:libfftw3-3.def" 
    & lib "/machine:X64" "/def:libfftw3f-3.def" 
    & lib "/machine:X64" "/def:libfftw3l-3.def" 
}
else{
    Write-Host "fftw seems to be installed already." 
}
Write-Host "" 
Write-Host "" 

cd $currentFolder 

#region Boost

if ((Test-Path -path "C:\PROGRA~1\boost\boost_1_67_0") -ne $True){   

    Write-Warning "Boost is about to be installed and compiled. This will take at least one hour... Take a break!" 
    Write-Warning "The script will wait after the Boost installation is completed.                               " 
    Write-Host -NoNewLine 'Press any key to continue...';
    Write-Host "" 
    Write-Host "" 

    $null = $Host.UI.RawUI.ReadKey('NoEcho,IncludeKeyDown');

    $boostFolder = $download_location + "\boost_binaries\boost_1_67_0" 

    #Download and install boost
    Write-Host "Downloading BOOOST (1.67.0)..." 
    download_file "https://dl.bintray.com/boostorg/release/1.67.0/source/boost_1_67_0.zip" ($download_location + "\boost_1_67_0.zip")
    Write-Host "Unzipping BOOST..." 
    unzip ($download_location + "\boost_1_67_0.zip")  "$download_location\boost_binaries" 

    cd $boostFolder

    Write-Host  "Building BJAM Binaries..." 

    if (Test-Path ".\b2.exe" -PathType Leaf){ del ".\b2.exe"}
    if (Test-Path ".\bjam.exe" -PathType Leaf){ del ".\bjam.exe"}

    & .\bootstrap.bat

    if (Test-Path ".\bin.v2"){ 
        rmdir .\bin.v2 /s /q
    }
    if (Test-Path "..\bin"){
        rmdir ..\bin /s /q
    }

    Write-Host "Building BOOST..." 

    if (Test-Path ".\stage_x64"){ 
        rmdir .\stage_x64 /s /q
    }

    &.\b2.exe --toolset=msvc-14.1 --clean-all
    &.\b2.exe --toolset=msvc-14.1 architecture=x86 address-model=64 --stagedir=".\stage_x64" threading=multi --build-type=complete stage

    cd $currentFolder

    #Make folders
    if (-Not (Test-Path "C:\PROGRA~1\boost")){ 
        mkdir "C:\PROGRA~1\boost" 
    }
    if (-Not (Test-Path "C:\PROGRA~1\boost\boost_1_67_0")){ 
        mkdir "C:\PROGRA~1\boost\boost_1_67_0" 
    }

    #Move files to Program folder
    Move-Item ($boostFolder+"\boost*\") "C:\PROGRA~1\boost\boost_1_67_0\boost" 
    Move-Item ($boostFolder+"\stage_x64\lib*\") "C:\PROGRA~1\boost\boost_1_67_0\lib" 
    Write-Host "" 
    Write-Host "" 

}
else{
    Write-Host "Boost (1.67.0) seems to be installed already.  " 
    Write-Host "" 
    Write-Host "" 
}
#endregion



###################################################################################
### Update system variables
###################################################################################
#region 
foreach($level in "Machine","User") {
   [Environment]::GetEnvironmentVariables($level).GetEnumerator() | % {
      # For Path variables, append the new values, if they're not already in there
      if($_.Name -match 'Path$') { 
         $_.Value = ($((Get-Content "Env:$($_.Name)") + ";$($_.Value)") -split ';' | Select -unique) -join ';'
      }
      $_
   } | Set-Content -Path { "Env:$($_.Name)" }
}
#endregion    

#region PATHS

    #Message reminding to set paths
    Write-Host "Ensuring that paths to the following locations are in the PATH environment variable: " 
    Write-Host "    - Boost libraries    (typically C:\Program Files\boost\boost_1_67\lib)" 
    Write-Host "    - FFTW libraries     (typically C:\MRILibraries\fftw3)" 
    Write-Host "    - HDF5 libraries     (typically C:\Program Files\HDF Group\HDF5\1.8.9\bin)" 
    Write-Host "" 
    Write-Host "" 

    ### Modify user environment variable
    [Environment]::GetEnvironmentVariables("User")     
    $pathToAdd = "$library_location\fftw3" 
    if ( -Not((Get-ItemProperty HKCU:\Environment).PATH  -like "*$pathToAdd3*")) {
        [Environment]::SetEnvironmentVariable("Path", (Get-ItemProperty HKCU:\Environment).PATH + ";$pathToAdd", [System.EnvironmentVariableTarget]::User)
    }
    $pathToAdd = "C:\Program Files\boost\boost_1_67\lib" 
    if ( -Not((Get-ItemProperty HKCU:\Environment).PATH  -like "*$pathToAdd*")) {
        [Environment]::SetEnvironmentVariable("Path", (Get-ItemProperty HKCU:\Environment).PATH + ";$pathToAdd", [System.EnvironmentVariableTarget]::User)
    }

    $pathToAdd = "C:\Program Files\HDF Group\HDF5\1.8.9\bin" 
    if ( -Not((Get-ItemProperty HKCU:\Environment).PATH  -like "*$pathToAdd*")) {
        [Environment]::SetEnvironmentVariable("Path", (Get-ItemProperty HKCU:\Environment).PATH + ";$pathToAdd", [System.EnvironmentVariableTarget]::User)
    }
    Write-Host "" 
    Write-Host "" 

#endregion    

#region Update system variables
foreach($level in "Machine","User") {
   [Environment]::GetEnvironmentVariables($level).GetEnumerator() | % {
      # For Path variables, append the new values, if they're not already in there
      if($_.Name -match 'Path$') { 
         $_.Value = ($((Get-Content "Env:$($_.Name)") + ";$($_.Value)") -split ';' | Select -unique) -join ';'
      }
      $_
   } | Set-Content -Path { "Env:$($_.Name)" }
}
#endregion    



###################################################################################
### Install the ISMRMRD format
###################################################################################
#region 

    Write-Host 'Now we are ready to install the ISMRMRD Format. Press any key to continue...';
    $null = $Host.UI.RawUI.ReadKey('NoEcho,IncludeKeyDown');
    Write-Host "" 
    Write-Host "" 

    #Now check if ISMRMRD folder already exists
    if ((Test-Path -path $library_location\ismrmrd)){ 
        Write-Warning "Folder $library_location\ismrmrd already exists. It will be removed.";
        $confirmation = Read-Host "Are you sure you want to proceed? [Y] Yes  [N] No (default is 'Y')" 
        if ($confirmation -eq 'N') {
          Write-Warning "Installation has been cancelled." 
            cd $currentFolder 
            break            
        }
        Remove-Item $library_location\ismrmrd -Force -Recurse
    }

    #Clone ISMRMRD repository
    Write-Host "Cloning ISMRMRD (master)..." 

    cd $library_location

    git clone "https://github.com/SkopeMagneticResonanceTechnologies/ismrmrd" 
	
    #Create solution
    if (-Not (Test-Path "$library_location\ismrmrd\build")){ 
        mkdir "$library_location\ismrmrd\build" 
    }
    cd "$library_location\ismrmrd\build" 
    Write-Host "Creating solution..." 
    &cmake ..\ "-G" "Visual Studio 15 2017 Win64" 

    #Cross your fingers. Here we go!
    Write-Host "Building ISMRMRD..." 
    &msbuild .\ISMRMRD.sln /p:Configuration='Release'

    Copy-Item "$library_location\ismrmrd\build\include\ismrmrd\version.h" "$library_location\ismrmrd\include\ismrmrd" 

    cd $currentFolder 

    $pathToAdd = "$library_location\ismrmrd\build\Release\" 
    if ( -Not((Get-ItemProperty HKCU:\Environment).PATH  -like "*$pathToAdd*")) {
        [Environment]::SetEnvironmentVariable("Path", (Get-ItemProperty HKCU:\Environment).PATH + ";$pathToAdd", [System.EnvironmentVariableTarget]::User)
    }

    Write-Host "" 
    Write-Host "" 

#endregion



###################################################################################
### Install dependencies
###################################################################################
#region XML, XSL stuff

    function download_unzip_file($app,$version) {  
	
		# Clean first
		if ((Test-Path -path ($download_location + "\$app-$version.7z"))){ 
			Write-Host "Folder $library_location\siemens_to_ismrmrd already exists. It will be removed.";
			Remove-Item "$library_location\$app-$version.7z" -Force -Recurse
		}
		if ((Test-Path -path $library_location\$app)){ 
			Write-Host "Folder $library_location\siemens_to_ismrmrd already exists. It will be removed.";
			Remove-Item $library_location\$app -Force -Recurse
		}
		
		# Download
		Write-Host "Downloading $app (iconv-$version)..."
        download_file "ftp://ftp.zlatkovic.com/pub/libxml/64bit/$app-$version.7z" ($download_location + "\$app-$version.7z")
		# Unzip
		
        Write-Host "Unzipping $app..." 
        sz x "$download_location\$app-$version.7z"  -o"$download_location\$app" 
		
		# Move 
        Move-Item -Path "$download_location\$app"  -Destination "$library_location\" -force
		
		# Add to path
        $pathToAdd = "$library_location\$app\bin" 
        if ( -Not((Get-ItemProperty HKCU:\Environment).PATH  -like "*$pathToAdd*")) {
            [Environment]::SetEnvironmentVariable("Path", (Get-ItemProperty HKCU:\Environment).PATH + ";$pathToAdd", [System.EnvironmentVariableTarget]::User)
        }
        Write-Host "" 
        Write-Host "" 
    }

    #Download unzip, and install iconv
    $app = "iconv" 
    $version = "1.14-win32-x86_64" 
    download_unzip_file "$app" "$version" 

    #Download, unzip, and install libxml2
    $app = "libxml2" 
    $libraryName = "libxml2-2" 
    $version = "2.9.3-win32-x86_64" 
    download_unzip_file "$app" "$version" 

    #Create lib file
    &Get-Exports -DllPath "$library_location\$app\bin\$libraryName.dll" -ExportsToCpp "$library_location\$app\bin\$libraryName.def" 
    Write-Host "Creating lib file: $library_location\$app\bin\$libraryName.dll" 
    Write-Host "" 
    Write-Host "" 
    &lib /def:"$library_location\$app\bin\$libraryName.def" /out:"$library_location\$app\lib\$libraryName.lib" 

    #Download, unzip, and install libxslt
    $app = "libxslt" 
    $libraryName = "libxslt-1" 
    $version = "1.1.28-win32-x86_64" 
    download_unzip_file "$app" "$version" 

    #Create lib file
    &Get-Exports -DllPath "$library_location\$app\bin\$libraryName.dll" -ExportsToCpp "$library_location\$app\bin\$libraryName.def" 
    Write-Host "Creating lib file: $library_location\$app\bin\$libraryName.dll" 
    Write-Host "" 
    Write-Host "" 
    &lib /def:"$library_location\$app\bin\$libraryName.def" /out:"$library_location\$app\lib\$libraryName.lib" 

    #Download, unzip, and install mingwrt
    $app = "mingwrt" 
    $version = "5.2.0-win32-x86_64" 
    download_unzip_file "$app" "$version"     

    #Download, unzip, and install zlib
    $app = "zlib" 
    $version = "1.2.8-win32-x86_64" 
    download_unzip_file "$app" "$version"  

	Write-Host "" 
    Write-Host ""	

#endregion



###################################################################################
### Update system variables
###################################################################################
#region
foreach($level in "Machine","User") {
   [Environment]::GetEnvironmentVariables($level)
}

foreach($level in "Machine","User") {
   [Environment]::GetEnvironmentVariables($level).GetEnumerator() | % {
      # For Path variables, append the new values, if they're not already in there
      if($_.Name -match 'Path$') { 
         $_.Value = ($((Get-Content "Env:$($_.Name)") + ";$($_.Value)") -split ';' | Select -unique) -join ';'
      }
      $_
   } | Set-Content -Path { "Env:$($_.Name)" }
}
#endregion    



###################################################################################
### Install the SIEMENS_TO_ISMRMRD converter
###################################################################################
#region 

    Write-Host 'Now we are ready to install the Siemens converter. Press any key to continue...';
    $null = $Host.UI.RawUI.ReadKey('NoEcho,IncludeKeyDown');
    Write-Host "" 
    Write-Host "" 

    #Now check if siemens_to_ismrmrd folder already exists
    if ((Test-Path -path $library_location\siemens_to_ismrmrd)){ 
        Write-Warning "Folder $library_location\siemens_to_ismrmrd already exists. It will be removed.";
        $confirmation = Read-Host "Are you sure you want to proceed? [Y] Yes  [N] No (default is 'Y')" 
        if ($confirmation -eq 'N') {
          Write-Warning "Installation has been cancelled." 
            cd $currentFolder 
            break            
        }
        Remove-Item $library_location\siemens_to_ismrmrd -Force -Recurse
    }

    #Now download SIEMENS_TO_ISMRMRD
    Write-Host "Cloning SIEMENS_TO_ISMRMRD (master)..." 
    cd $library_location    
    git clone "https://github.com/SkopeMagneticResonanceTechnologies/siemens_to_ismrmrd" 

    cd "$library_location\siemens_to_ismrmrd" 
    #git reset --hard "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" # Reset to previous version

    #Create solution
    Write-Host "Creating solution..." 
    if (-Not (Test-Path "$library_location\siemens_to_ismrmrd\build")){ 
        mkdir "$library_location\siemens_to_ismrmrd\build" 
    }
    cd "$library_location\siemens_to_ismrmrd\build" 
    & cmake ..\ "-G" "Visual Studio 15 2017 Win64" "-DISMRMRD_DIR=$library_location\ismrmrd\build\InstallFiles" "-DLIBXML2_LIBRARY=$library_location\libxml2\lib\libxml2-2.lib" "-DLIBXSLT_LIBRARIES=$library_location\libxslt\lib\libxslt-1.lib" "-DLIBXSLT_LIBRARIES=$library_location\libxslt\lib\libxslt-1.lib" "-DCMAKE_CXX_FLAGS=-I$library_location\iconv\include" 

    #Cross your fingers. Here we go!
    Write-Host "Building SIEMENS_TO_ISMRMRD..." 
    &msbuild .\SIEMENS-TO-ISMRMRD.sln /p:Configuration='Release' 

    cd $currentFolder 

    Write-Host "" 
    Write-Host "" 

    Write-Host "Copying binaries to C:\PROGRA~1\SIEMENS_TO_ISMRMRD..." 
    Write-Host "" 
    Write-Host "" 
    if (-Not (Test-Path "C:\PROGRA~1\SIEMENS_TO_ISMRMRD")){ 
        mkdir "C:\PROGRA~1\SIEMENS_TO_ISMRMRD" 
    }

    $siemensConverter_location = "C:\PROGRA~1\SIEMENS_TO_ISMRMRD\" 

    Copy-Item -Path "$library_location\siemens_to_ismrmrd\build\Release\siemens_to_ismrmrd.exe" -Destination $siemensConverter_location  -force
    Copy-Item -Path "$library_location\siemens_to_ismrmrd\build\Release\embed.exe" -Destination $siemensConverter_location -force

    Write-Host "Adding other dlls to C:\PROGRA~1\SIEMENS_TO_ISMRMRD..." 
    Write-Host "" 
    Write-Host "" 
    Copy-Item -Path "$library_location\iconv\bin\libiconv-2.dll" -Destination $siemensConverter_location -force    
    Copy-Item -Path "$library_location\libxml2\bin\libxml2-2.dll" -Destination $siemensConverter_location -force    
    Copy-Item -Path "$library_location\libxslt\bin\libxslt-1.dll" -Destination $siemensConverter_location -force    
    Copy-Item -Path "$library_location\libxslt\bin\xsltproc.exe" -Destination $siemensConverter_location -force    
    Copy-Item -Path "$library_location\mingwrt\bin\libwinpthread-1.dll" -Destination $siemensConverter_location -force    
    Copy-Item -Path "$library_location\zlib\bin\zlib1.dll" -Destination $siemensConverter_location -force    

    Write-Host "Adding ismrmrd.dll to C:\PROGRA~1\SIEMENS_TO_ISMRMRD..." 
    Write-Host "" 
    Write-Host "" 
    Copy-Item -Path "$library_location\ismrmrd\build\Release\ismrmrd.dll" -Destination $siemensConverter_location -force

    Write-Host "Adding SIEMENS_TO_ISMRMRD to PATH..." 
    Write-Host "" 
    Write-Host "" 
    $pathToAdd = $siemensConverter_location
    if ( -Not((Get-ItemProperty HKCU:\Environment).PATH  -like "*$pathToAdd*")) {
        [Environment]::SetEnvironmentVariable("Path", (Get-ItemProperty HKCU:\Environment).PATH + ";$pathToAdd", [System.EnvironmentVariableTarget]::User)
    }

#endregion



###################################################################################
### Cleanup
###################################################################################
Remove-Item $download_location -Force -Recurse
Write-Host "**********************************************************************" 
Write-Host "**********************************************************************" 
Write-Host "** HOPEFULLY EVERTHING WORKS NOW                                      " 
Write-Host "**********************************************************************" 
Write-Host "**********************************************************************" 