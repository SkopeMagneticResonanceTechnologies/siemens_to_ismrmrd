<?xml version="1.0" encoding="ISO-8859-1"?>

<xsl:stylesheet version="1.0"
xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

    <xsl:output method="xml" indent="yes"/>

    <xsl:variable name="phaseOversampling">
        <xsl:choose>
            <xsl:when test="not(siemens/IRIS/DERIVED/phaseOversampling)">0</xsl:when>
            <xsl:otherwise>
                <xsl:value-of select="siemens/IRIS/DERIVED/phaseOversampling"/>
            </xsl:otherwise>
        </xsl:choose>
    </xsl:variable>

    <xsl:variable name="sliceOversampling">
        <xsl:choose>
            <xsl:when test="not(siemens/MEAS/sKSpace/dSliceOversamplingForDialog)">0</xsl:when>
            <xsl:otherwise>
                <xsl:value-of select="siemens/MEAS/sKSpace/dSliceOversamplingForDialog"/>
            </xsl:otherwise>
        </xsl:choose>
    </xsl:variable>

    <xsl:variable name="partialFourierPhase">
        <xsl:choose>
            <xsl:when test="siemens/MEAS/sKSpace/ucPhasePartialFourier = 1">0.5</xsl:when>
            <xsl:when test="siemens/MEAS/sKSpace/ucPhasePartialFourier = 2">0.625</xsl:when>
            <xsl:when test="siemens/MEAS/sKSpace/ucPhasePartialFourier = 4">0.75</xsl:when>
            <xsl:when test="siemens/MEAS/sKSpace/ucPhasePartialFourier = 8">0.875</xsl:when>
            <xsl:otherwise>1.0</xsl:otherwise>
        </xsl:choose>
    </xsl:variable>
	
	<xsl:variable name="readoutOversampling">
		<xsl:choose>
			<xsl:when test="not(siemens/YAPS/flReadoutOSFactor)">1</xsl:when>
			<xsl:otherwise>
				<xsl:value-of select="siemens/YAPS/flReadoutOSFactor"/>
			</xsl:otherwise>
		</xsl:choose>
	</xsl:variable>
	
	<xsl:variable name="accelerationFactor1">
		<xsl:choose>
			<xsl:when test="not(siemens/MEAS/sPat/lAccelFactPE)">1</xsl:when>
			<xsl:otherwise>
				<xsl:value-of select="(siemens/MEAS/sPat/lAccelFactPE)"/>
			</xsl:otherwise>
		</xsl:choose>
	</xsl:variable>

	<xsl:variable name="accelerationFactor2">
		<xsl:choose>
			<xsl:when test="not(siemens/MEAS/sPat/lAccelFact3D)">1</xsl:when>
			<xsl:otherwise>
				<xsl:value-of select="(siemens/MEAS/sPat/lAccelFact3D)"/>
			</xsl:otherwise>
		</xsl:choose>
	</xsl:variable>
	

    <xsl:variable name="numberOfContrasts">
        <xsl:value-of select="siemens/MEAS/lContrasts"/>
    </xsl:variable>

    <xsl:variable name="studyID">
        <xsl:value-of select="substring(siemens/IRIS/RECOMPOSE/StudyLOID, 6)"/>
    </xsl:variable>

    <xsl:variable name="patientID">
        <xsl:value-of select="substring(siemens/IRIS/RECOMPOSE/PatientLOID, 6)"/>
    </xsl:variable>
		
	
	<xsl:variable name="birthDate">	
		<xsl:choose>
		  <xsl:when test="siemens/IRIS/RECOMPOSE/PatientBirthDay='xxxxxxxx'">
			<xsl:value-of select="19000101"/>
		  </xsl:when>
		  <xsl:otherwise>
			<xsl:value-of select="siemens/IRIS/RECOMPOSE/PatientBirthDay"/>
		  </xsl:otherwise>
		</xsl:choose> 
	</xsl:variable>
		
	<xsl:variable name="birthYear">
        <xsl:value-of select="substring($birthDate,1,4)"/>
    </xsl:variable>
	<xsl:variable name="birthMonth">
        <xsl:value-of select="substring($birthDate,5,2)"/>
    </xsl:variable>
	<xsl:variable name="birthDay">
        <xsl:value-of select="substring($birthDate,7,2)"/>
    </xsl:variable>
	
	<!-- The date of the study is given after the 10th dot in the tFrameOfReference -->
	<xsl:variable name="temp1"><xsl:value-of select="substring-after(string(siemens/YAPS/tFrameOfReference), '.')" /></xsl:variable> 
	<xsl:variable name="temp2"><xsl:value-of select="substring-after($temp1, '.')" /></xsl:variable> 
	<xsl:variable name="temp3"><xsl:value-of select="substring-after($temp2, '.')" /></xsl:variable> 
	<xsl:variable name="temp4"><xsl:value-of select="substring-after($temp3, '.')" /></xsl:variable> 
	<xsl:variable name="temp5"><xsl:value-of select="substring-after($temp4, '.')" /></xsl:variable> 
	<xsl:variable name="temp6"><xsl:value-of select="substring-after($temp5, '.')" /></xsl:variable> 
	<xsl:variable name="temp7"><xsl:value-of select="substring-after($temp6, '.')" /></xsl:variable> 
	<xsl:variable name="temp8"><xsl:value-of select="substring-after($temp7, '.')" /></xsl:variable> 
	<xsl:variable name="temp9"><xsl:value-of select="substring-after($temp8, '.')" /></xsl:variable> 
	<xsl:variable name="studyDate"><xsl:value-of select="substring-after($temp9, '.')" /></xsl:variable> 
	 
	<xsl:variable name="studyYear">
        <xsl:value-of select="substring($studyDate,1,4)"/>
    </xsl:variable>
	<xsl:variable name="studyMonth">
        <xsl:value-of select="substring($studyDate,5,2)"/>
    </xsl:variable>
	<xsl:variable name="studyDay">
        <xsl:value-of select="substring($studyDate,7,2)"/>
    </xsl:variable>
	
	<xsl:variable name="dateSeperator">-</xsl:variable>    
    <xsl:variable name="strSeperator">_</xsl:variable>

    <xsl:template match="/">
        <ismrmrdHeader xsi:schemaLocation="http://www.ismrm.org/ISMRMRD ismrmrd.xsd"
                xmlns="http://www.ismrm.org/ISMRMRD"
                xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
                xmlns:xs="http://www.w3.org/2001/XMLSchema">

			<subjectInformation>
                <patientName>
                    <xsl:value-of select="siemens/DICOM/tPatientName"/>
                </patientName>
                <xsl:if test="siemens/YAPS/flUsedPatientWeight > 0">
                    <patientWeight_kg>
                        <xsl:value-of select="siemens/YAPS/flUsedPatientWeight"/>
                    </patientWeight_kg>
                </xsl:if>
								
				<patientID>
					<xsl:value-of select="siemens/IRIS/RECOMPOSE/PatientID"/>
				</patientID> 
				<patientBirthdate>
					<xsl:value-of select="concat($birthYear,$dateSeperator,$birthMonth,$dateSeperator,$birthDay)"/>
				</patientBirthdate> 
                <patientGender>
                    <xsl:choose>
                        <xsl:when test="siemens/DICOM/lPatientSex = 1">F</xsl:when>
                        <xsl:when test="siemens/DICOM/lPatientSex = 2">M</xsl:when>
                        <xsl:otherwise>O</xsl:otherwise>
                    </xsl:choose>
                </patientGender>
            </subjectInformation>
			
			<studyInformation>
				<studyInstanceUID>
                    <xsl:value-of select="$studyID" />
                </studyInstanceUID>
				<studyDate>
					<xsl:value-of select="concat($studyYear,$dateSeperator,$studyMonth,$dateSeperator,$studyDay)"/>
				</studyDate>
            </studyInformation>

            <measurementInformation>
                <measurementID>
                    <xsl:value-of select="concat(string(siemens/DICOM/DeviceSerialNumber), $strSeperator, $patientID, $strSeperator, $studyID, $strSeperator, string(siemens/HEADER/MeasUID))"/>
                </measurementID>
                <patientPosition>
                    <xsl:value-of select="siemens/YAPS/tPatientPosition"/>
                </patientPosition>
                <protocolName>
                    <xsl:value-of select="siemens/MEAS/tProtocolName"/>
                </protocolName>

                <xsl:if test="siemens/YAPS/ReconMeasDependencies/RFMap > 0">
                    <measurementDependency>
                        <dependencyType>RFMap</dependencyType>
                        <measurementID>
                            <xsl:value-of select="concat(string(siemens/DICOM/DeviceSerialNumber), $strSeperator, $patientID, $strSeperator, $studyID, $strSeperator, string(siemens/YAPS/ReconMeasDependencies/RFMap))"/>
                        </measurementID>
                    </measurementDependency>
                </xsl:if>

                <xsl:if test="siemens/YAPS/ReconMeasDependencies/SenMap > 0">
                    <measurementDependency>
                        <dependencyType>SenMap</dependencyType>
                        <measurementID>
                            <xsl:value-of select="concat(string(siemens/DICOM/DeviceSerialNumber), $strSeperator, $patientID, $strSeperator, $studyID, $strSeperator, string(siemens/YAPS/ReconMeasDependencies/SenMap))"/>
                        </measurementID>
                    </measurementDependency>
                </xsl:if>

                <xsl:if test="siemens/YAPS/ReconMeasDependencies/Noise > 0">
                    <measurementDependency>
                        <dependencyType>Noise</dependencyType>
                        <measurementID>
                            <xsl:value-of select="concat(string(siemens/DICOM/DeviceSerialNumber), $strSeperator, $patientID, $strSeperator, $studyID, $strSeperator, string(siemens/YAPS/ReconMeasDependencies/Noise))"/>
                        </measurementID>
                    </measurementDependency>
                </xsl:if>

                <frameOfReferenceUID>
                    <xsl:value-of select="siemens/YAPS/tFrameOfReference" />
                </frameOfReferenceUID>

            </measurementInformation>

            <acquisitionSystemInformation>
                <systemVendor>
                    <xsl:value-of select="siemens/DICOM/Manufacturer"/>
                </systemVendor>
                <systemModel>
                    <xsl:value-of select="siemens/DICOM/ManufacturersModelName"/>
                </systemModel>
                <systemFieldStrength_T>
                    <xsl:value-of select="siemens/YAPS/flMagneticFieldStrength"/>
                </systemFieldStrength_T>
                <relativeReceiverNoiseBandwidth>0.793</relativeReceiverNoiseBandwidth>
                <receiverChannels>
                    <xsl:value-of select="siemens/YAPS/iMaxNoOfRxChannels" />
                </receiverChannels>

		<!-- Coil Labels -->
		<xsl:choose>
                  <!-- VD line with dual density -->
                    <xsl:when test="siemens/MEAS/asCoilSelectMeas/ADC/lADCChannelConnected">
                        <xsl:variable name="NumberOfSelectedCoils">
                            <xsl:value-of select="count(siemens/MEAS/asCoilSelectMeas/Select/lElementSelected[text() = '1'])" />
                        </xsl:variable>
                        <xsl:for-each select="siemens/MEAS/asCoilSelectMeas/ADC/lADCChannelConnected[position() >= 1  and not(position() > $NumberOfSelectedCoils)]">
                            <xsl:sort data-type="number"
                                      select="." />
                            <xsl:variable name="CurADC"
                                          select="."/>
                            <xsl:variable name="CurADCIndex"
                                          select="position()" />
                            <xsl:for-each select="../lADCChannelConnected[position() >= 1  and not(position() > $NumberOfSelectedCoils)]">
                                <xsl:if test="$CurADC = .">
                                    <xsl:variable name="CurCoil" select="position()"/>
                                    <xsl:variable name="CurCoilID" select="../../ID/tCoilID[$CurCoil]"/>
                                    <xsl:variable name="CurCoilElement" select="../../Elem/tElement[$CurCoil]"/>
                                    <xsl:variable name="CurCoilCopyID" select="../../Coil/lCoilCopy[$CurCoil]"/>
                                    <coilLabel>
                                        <coilNumber>
                                            <xsl:value-of select="number(../lADCChannelConnected[$CurADCIndex])"/>
                                        </coilNumber>
                                        <coilName>
                                            <xsl:value-of select="$CurCoilID"/>:<xsl:value-of select="string($CurCoilCopyID)"/>:<xsl:value-of select="$CurCoilElement"/>
                                        </coilName>
                                    </coilLabel>
                                </xsl:if>
                            </xsl:for-each>
                        </xsl:for-each>
                    </xsl:when>
                  <xsl:otherwise> <!-- This is probably VB -->
                    <xsl:for-each select="siemens/MEAS/asCoilSelectMeas/ID/tCoilID">
                      <xsl:variable name="CurCoil" select="position()"/>
		      <coilLabel>
			<coilNumber><xsl:value-of select="$CurCoil -1"/></coilNumber>
			<coilName><xsl:value-of select="."/>:<xsl:value-of select="../../Elem/tElement[$CurCoil]"/></coilName>
		      </coilLabel>
                    </xsl:for-each>
                  </xsl:otherwise>
                </xsl:choose>

                <institutionName>
                    <xsl:value-of select="siemens/DICOM/InstitutionName" />
                </institutionName>
            </acquisitionSystemInformation>

            <experimentalConditions>
                <H1resonanceFrequency_Hz>
                    <xsl:value-of select="siemens/DICOM/lFrequency"/>
                </H1resonanceFrequency_Hz>
				<LarmorConstant_Hz_Per_T>
                    <xsl:value-of select="siemens/YAPS/alLarmorConstant"/>
                </LarmorConstant_Hz_Per_T>
            </experimentalConditions>

            <encoding>
                <trajectory>spiral</trajectory>
                <trajectoryDescription>
                    <identifier>spiral</identifier>				
				
					
                    <userParameterDouble>
                        <name>dwellTime</name>
                        <value>
                            <xsl:value-of select="siemens/MEAS/sRXSPEC/alDwellTime div 1000.0"/>
                        </value>
                    </userParameterDouble>
										
                    <comment>Spiral sequence</comment>
                </trajectoryDescription>
				
				<!--Hardcoded in sequence  -->
				<echoTrainLength>10</echoTrainLength>
				
                <encodedSpace>
                    <matrixSize>
						<!--Hardcoded value in sequence-->
						<x>2000</x>
						<y>
							<xsl:value-of select="floor(siemens/YAPS/iNoOfFourierLines div siemens/YAPS/GradientEchoTrainLength)"/>
						</y>
						
						<xsl:choose>
							<xsl:when test="not(siemens/YAPS/iNoOfFourierPartitions) or (siemens/YAPS/i3DFTLength = 1)">
								<z>1</z>
							</xsl:when>
							<xsl:otherwise>
								<z>
									<xsl:value-of select="siemens/YAPS/i3DFTLength"/>
								</z>
							</xsl:otherwise>
						</xsl:choose>
                    </matrixSize>
                    <fieldOfView_mm>
                        <x>
                            <xsl:value-of select="siemens/MEAS/sSliceArray/asSlice/s0/dReadoutFOV * $readoutOversampling"/>
                        </x>
                        <y>
                            <xsl:value-of select="siemens/MEAS/sSliceArray/asSlice/s0/dPhaseFOV * (1+$phaseOversampling)"/>
                        </y>
                        <z>
                            <xsl:value-of select="siemens/MEAS/sSliceArray/asSlice/s0/dThickness * (1+$sliceOversampling)"/>
                        </z>
                    </fieldOfView_mm>
                </encodedSpace>
                <reconSpace>
                    <matrixSize>
                        <x>
                            <xsl:value-of select="siemens/IRIS/DERIVED/imageColumns"/>
                        </x>
                        <y>
                            <xsl:value-of select="siemens/IRIS/DERIVED/imageLines"/>
                        </y>
                        <xsl:choose>
                            <xsl:when test="siemens/YAPS/i3DFTLength = 1">
                                <z>1</z>
                            </xsl:when>
                            <xsl:otherwise>
                                <z>
                                    <xsl:value-of select="siemens/MEAS/sKSpace/lImagesPerSlab"/>
                                </z>
                            </xsl:otherwise>
                        </xsl:choose>
                    </matrixSize>
                    <fieldOfView_mm>
                        <x>
                            <xsl:value-of select="siemens/MEAS/sSliceArray/asSlice/s0/dReadoutFOV"/>
                        </x>
                        <y>
                            <xsl:value-of select="siemens/MEAS/sSliceArray/asSlice/s0/dPhaseFOV"/>
                        </y>
                        <z>
                            <xsl:value-of select="siemens/MEAS/sSliceArray/asSlice/s0/dThickness"/>
                        </z>
                    </fieldOfView_mm>
                </reconSpace>
                <encodingLimits>
                    <kspace_encoding_step_1>
                        <minimum>0</minimum>
                        <maximum>
                            <xsl:value-of select="floor(siemens/YAPS/iNoOfFourierLines div siemens/YAPS/GradientEchoTrainLength) - 1"/>
                        </maximum>
                        <center>0</center>
                    </kspace_encoding_step_1>
                    <kspace_encoding_step_2>
                        <minimum>0</minimum>
						<maximum>0</maximum>
						<center>0</center>
                    </kspace_encoding_step_2>
                    <slice>
                        <minimum>0</minimum>
                        <maximum>
                            <xsl:value-of select="siemens/MEAS/sSliceArray/lSize - 1"/>
                        </maximum>
                        <center>0</center>
                    </slice>
                    <set>
                        <minimum>0</minimum>
                        <maximum>
                            <xsl:choose>
                                <xsl:when test="siemens/MEAS/sDiffusion/lDiffWeightings">
                                    <xsl:value-of select="siemens/MEAS/sDiffusion/lDiffWeightings - 1"/>
                                </xsl:when>
                                <xsl:otherwise>0</xsl:otherwise>
                            </xsl:choose>
                        </maximum>
                        <center>0</center>
                    </set>
                    <phase>
                        <minimum>0</minimum>
                        <maximum>
                            <xsl:choose>
                                <xsl:when test="siemens/MEAS/sPhysioImaging/lPhases">
                                    <xsl:value-of select="siemens/MEAS/sPhysioImaging/lPhases - 1"/>
                                </xsl:when>
                                <xsl:otherwise>0</xsl:otherwise>
                            </xsl:choose>
                        </maximum>
                        <center>0</center>
                    </phase>
                    <repetition>
                        <minimum>0</minimum>
                        <maximum>
                            <xsl:choose>
                                <xsl:when test="siemens/MEAS/lRepetitions">
                                    <xsl:value-of select="siemens/MEAS/lRepetitions"/>
                                </xsl:when>
                                <xsl:otherwise>0</xsl:otherwise>
                            </xsl:choose>
                        </maximum>
                        <center>0</center>
                    </repetition>
                    <segment>
                        <minimum>0</minimum>
						<!--Hardcoded value in sequence-->
                        <maximum>9</maximum>
                        <center>0</center>
                    </segment>
                    <contrast>
                        <minimum>0</minimum>
                        <maximum>
                            <xsl:choose>
                                <xsl:when test="siemens/MEAS/lContrasts">
                                    <xsl:value-of select="siemens/MEAS/lContrasts - 1"/>
                                </xsl:when>
                                <xsl:otherwise>0</xsl:otherwise>
                            </xsl:choose>
                        </maximum>
                        <center>0</center>
                    </contrast>
                    <average>
                        <minimum>0</minimum>
                        <maximum>
                            <xsl:choose>
                                <xsl:when test="siemens/MEAS/lAverages">
                                    <xsl:value-of select="siemens/MEAS/lAverages - 1"/>
                                </xsl:when>
                                <xsl:otherwise>0</xsl:otherwise>
                            </xsl:choose>
                        </maximum>
                        <center>0</center>
                    </average>
					<user_0>
						<minimum>0</minimum>
                        <maximum>
                            <xsl:choose>
                                <xsl:when test="siemens/MEAS/sDiffusion/lDiffDirections">
                                    <xsl:value-of select="siemens/MEAS/sDiffusion/lDiffDirections - 1"/>
                                </xsl:when>
                                <xsl:otherwise>0</xsl:otherwise>
                            </xsl:choose>
                        </maximum>
                        <center>0</center>
					</user_0>
					<user_1>
						<minimum>0</minimum>
						<maximum>0</maximum>
						<center>0</center>
					</user_1>
					<user_2>
						<minimum>0</minimum>
						<maximum>0</maximum>
						<center>0</center>
					</user_2>
					<user_3>
						<minimum>0</minimum>
						<maximum>0</maximum>
						<center>0</center>
					</user_3>
					<user_4>
						<minimum>0</minimum>
						<maximum>0</maximum>
						<center>0</center>
					</user_4>
					<user_5>
						<minimum>0</minimum>
						<maximum>0</maximum>
						<center>0</center>
					</user_5>
					<user_6>
						<minimum>0</minimum>
						<maximum>0</maximum>
						<center>0</center>
					</user_6>
					<user_7>
						<minimum>0</minimum>
						<maximum>0</maximum>
						<center>0</center>
					</user_7>
                </encodingLimits>
                
                <parallelImaging>
                  <accelerationFactor>
                            <kspace_encoding_step_1>
                      <xsl:choose>
                    <xsl:when test="not(siemens/MEAS/sPat/lAccelFactPE)">1</xsl:when>
                    <xsl:otherwise>
                      <xsl:value-of select="(siemens/MEAS/sPat/lAccelFactPE)"/>
                    </xsl:otherwise>
                      </xsl:choose>
                            </kspace_encoding_step_1>
                            <kspace_encoding_step_2>
                      <xsl:choose>
                    <xsl:when test="not(siemens/MEAS/sPat/lAccelFact3D)">1</xsl:when>
                    <xsl:otherwise>
                      <xsl:value-of select="(siemens/MEAS/sPat/lAccelFact3D)"/>
                    </xsl:otherwise>
                      </xsl:choose>
                            </kspace_encoding_step_2>
                  </accelerationFactor>
                  <calibrationMode>
                            <xsl:choose>
                      <xsl:when test="siemens/MEAS/sPat/ucRefScanMode = 1">other</xsl:when>
                      <xsl:when test="siemens/MEAS/sPat/ucRefScanMode = 2">embedded</xsl:when>
                      <xsl:when test="siemens/MEAS/sPat/ucRefScanMode = 4">separate</xsl:when>
                      <xsl:when test="siemens/MEAS/sPat/ucRefScanMode = 8">separate</xsl:when>
                      <xsl:when test="siemens/MEAS/sPat/ucRefScanMode = 16">interleaved</xsl:when>
                      <xsl:when test="siemens/MEAS/sPat/ucRefScanMode = 32">interleaved</xsl:when>
                      <xsl:when test="siemens/MEAS/sPat/ucRefScanMode = 64">interleaved</xsl:when>
                      <xsl:otherwise>other</xsl:otherwise>
                            </xsl:choose>
                  </calibrationMode>
                  <xsl:if test="(siemens/MEAS/sPat/ucRefScanMode = 1) or (siemens/MEAS/sPat/ucRefScanMode = 16) or (siemens/MEAS/sPat/ucRefScanMode = 32) or (siemens/MEAS/sPat/ucRefScanMode = 64)">
                            <interleavingDimension>
                      <xsl:choose>
                    <xsl:when test="siemens/MEAS/sPat/ucRefScanMode = 16">average</xsl:when>
                    <xsl:when test="siemens/MEAS/sPat/ucRefScanMode = 32">repetition</xsl:when>
                    <xsl:when test="siemens/MEAS/sPat/ucRefScanMode = 64">phase</xsl:when>
                    <xsl:otherwise>other</xsl:otherwise>
                      </xsl:choose>
                            </interleavingDimension>
                  </xsl:if>
                </parallelImaging>
				
			
						
            </encoding>

            <sequenceParameters>
                <xsl:for-each select="siemens/MEAS/alTR">
                    <xsl:if test=". &gt; 0">
                        <TR>
                            <xsl:value-of select=". div 1000.0" />
                        </TR>
                    </xsl:if>
                </xsl:for-each>
                <xsl:for-each select="siemens/MEAS/alTE">
                    <xsl:if test=". &gt; 0">
                        <xsl:if test="position() &lt; ($numberOfContrasts + 1)">
                            <TE>
                                <xsl:value-of select=". div 1000.0" />
                            </TE>
                        </xsl:if>
                    </xsl:if>
                </xsl:for-each>
                <xsl:for-each select="siemens/MEAS/alTI">
                    <xsl:if test=". &gt; 0">
                        <TI>
                            <xsl:value-of select=". div 1000.0" />
                        </TI>
                    </xsl:if>
                </xsl:for-each>
                <xsl:for-each select="siemens/DICOM/adFlipAngleDegree">
                <xsl:if test=". &gt; 0">
                        <flipAngle_deg>
                            <xsl:value-of select="." />
                        </flipAngle_deg>
                    </xsl:if>
                </xsl:for-each>
				 <!--It is a spin-echo sequence-->
                <sequence_type>SE</sequence_type>
                <xsl:if test="siemens/YAPS/lEchoSpacing">
                    <echo_spacing>
                        <!--Hardcoded value in sequence and rounded to gradient raster time-->
                            <xsl:value-of select="ceiling(siemens/MEAS/sRXSPEC/alDwellTime div 1000.0 * 2000 div 10)*10 div 1000"/>
                    </echo_spacing>
                </xsl:if>
            </sequenceParameters>
			
			
			<userParameters>
					<!-- Decide whether the converter should undo field-of-view shifts. 
						 The undoing of FoVShift by the converter is suboptimal because of clock-shifts.
						 The FoV shift should be disabled on the scanner.
					     By default, the converter does not change the raw data -->
					<userParameterLong>
						<name>undoFieldOfViewShifts</name>
						<value>
							<xsl:value-of select="0"/>
						</value>
					</userParameterLong>
					
					<!--If the multi-slice mode is 2 (interleaved), the converter will remap the  
					slice indices to match the anatomical order -->
					<xsl:if test="siemens/MEAS/sKSpace/ucMultiSliceMode">
						<userParameterLong>
							<name>multiSliceMode</name>
							<value><xsl:value-of select="siemens/MEAS/sKSpace/ucMultiSliceMode"/></value>
						</userParameterLong>
					</xsl:if>
					
					<!-- Trig to ADC time: imaging scans -->
					<userParameterDouble>
					    <!--Hardcoded value in sequence-->
						<name>triggerToAcquisitionDelay</name>
						<value>0.200</value>
					</userParameterDouble>
								

					<!-- Suggested minimal TR for skope monitoring system -->
					<userParameterDouble>
						<name>skopeMinimalTR</name>
						<value>
							<xsl:value-of select="siemens/MEAS/sWipMemBlock/adFree[1]"/>
						</value>
					</userParameterDouble>

					<!-- Suggested interleave TR for skope monitoring system -->
					<userParameterDouble>
						<name>skopeInterleaveTR</name>
						<value>
							<xsl:value-of select="siemens/MEAS/sWipMemBlock/adFree[2]"/>
						</value>
					</userParameterDouble>

					<!-- Suggested skope field probe acquisition duration -->
					<userParameterDouble>
						<name>skopeAqDur</name>
						<!--Hardcoded value in sequence-->
						<value>32</value>
					</userParameterDouble>

					<!-- Gradient free interval for field probe excitation -->
					<userParameterDouble>
						<name>skopeGradientFreeInterval</name>
						<value>0.000200</value>
					</userParameterDouble>
					
					
					<!-- ECC  -->
					<!-- B0 compensation amplitude X -->
					<xsl:if test="siemens/MEAS/sGRADSPEC/sB0CompensationX/aflAmplitude[1]">
					  <userParameterDouble>
						  <name>sB0CompensationX_Amp_0</name>
						  <value>
							  <xsl:value-of select="siemens/MEAS/sGRADSPEC/sB0CompensationX/aflAmplitude[1]" />
						  </value>
					  </userParameterDouble>
					</xsl:if>
					
					<xsl:if test="siemens/MEAS/sGRADSPEC/sB0CompensationX/aflAmplitude[2]">
					  <userParameterDouble>
						  <name>sB0CompensationX_Amp_1</name>
						  <value>
							  <xsl:value-of select="siemens/MEAS/sGRADSPEC/sB0CompensationX/aflAmplitude[2]" />
						  </value>
					  </userParameterDouble>
					</xsl:if>
					
					<xsl:if test="siemens/MEAS/sGRADSPEC/sB0CompensationX/aflAmplitude[3]">
					  <userParameterDouble>
						  <name>sB0CompensationX_Amp_2</name>
						  <value>
							  <xsl:value-of select="siemens/MEAS/sGRADSPEC/sB0CompensationX/aflAmplitude[3]" />
						  </value>
					  </userParameterDouble>
					</xsl:if>
					
					<!-- B0 compensation amplitude Y -->
					<xsl:if test="siemens/MEAS/sGRADSPEC/sB0CompensationY/aflAmplitude[1]">
					  <userParameterDouble>
						  <name>sB0CompensationY_Amp_0</name>
						  <value>
							  <xsl:value-of select="siemens/MEAS/sGRADSPEC/sB0CompensationY/aflAmplitude[1]" />
						  </value>
					  </userParameterDouble>
					</xsl:if>
					
					
					<xsl:if test="siemens/MEAS/sGRADSPEC/sB0CompensationY/aflAmplitude[2]">
					  <userParameterDouble>
						  <name>sB0CompensationY_Amp_1</name>
						  <value>
							  <xsl:value-of select="siemens/MEAS/sGRADSPEC/sB0CompensationY/aflAmplitude[2]" />
						  </value>
					  </userParameterDouble>
					</xsl:if>
					
					
					<xsl:if test="siemens/MEAS/sGRADSPEC/sB0CompensationY/aflAmplitude[3]">
					  <userParameterDouble>
						  <name>sB0CompensationY_Amp_2</name>
						  <value>
							  <xsl:value-of select="siemens/MEAS/sGRADSPEC/sB0CompensationY/aflAmplitude[3]" />
						  </value>
					  </userParameterDouble>
					</xsl:if>
					
					
					<!-- B0 compensation amplitude Z -->
					<xsl:if test="siemens/MEAS/sGRADSPEC/sB0CompensationZ/aflAmplitude[1]">
					  <userParameterDouble>
						  <name>sB0CompensationZ_Amp_0</name>
						  <value>
							  <xsl:value-of select="siemens/MEAS/sGRADSPEC/sB0CompensationZ/aflAmplitude[1]" />
						  </value>
					  </userParameterDouble>
					</xsl:if>
					
					
					<xsl:if test="siemens/MEAS/sGRADSPEC/sB0CompensationZ/aflAmplitude[2]">
					  <userParameterDouble>
						  <name>sB0CompensationZ_Amp_1</name>
						  <value>
							  <xsl:value-of select="siemens/MEAS/sGRADSPEC/sB0CompensationZ/aflAmplitude[2]" />
						  </value>
					  </userParameterDouble>
					</xsl:if>
					
					
					<xsl:if test="siemens/MEAS/sGRADSPEC/sB0CompensationZ/aflAmplitude[3]">
					  <userParameterDouble>
						  <name>sB0CompensationZ_Amp_2</name>
						  <value>
							  <xsl:value-of select="siemens/MEAS/sGRADSPEC/sB0CompensationZ/aflAmplitude[3]" />
						  </value>
					  </userParameterDouble>
					</xsl:if>				
					
					<!-- B0 compensation constant X -->
					<xsl:if test="siemens/MEAS/sGRADSPEC/sB0CompensationX/aflTimeConstant[1]">
					  <userParameterDouble>
						  <name>sB0CompensationX_Tau_0</name>
						  <value>
							  <xsl:value-of select="siemens/MEAS/sGRADSPEC/sB0CompensationX/aflTimeConstant[1]" />
						  </value>
					  </userParameterDouble>
					</xsl:if>
					
					<xsl:if test="siemens/MEAS/sGRADSPEC/sB0CompensationX/aflTimeConstant[2]">
					  <userParameterDouble>
						  <name>sB0CompensationX_Tau_1</name>
						  <value>
							  <xsl:value-of select="siemens/MEAS/sGRADSPEC/sB0CompensationX/aflTimeConstant[2]" />
						  </value>
					  </userParameterDouble>
					</xsl:if>
					
					<xsl:if test="siemens/MEAS/sGRADSPEC/sB0CompensationX/aflTimeConstant[3]">
					  <userParameterDouble>
						  <name>sB0CompensationX_Tau_2</name>
						  <value>
							  <xsl:value-of select="siemens/MEAS/sGRADSPEC/sB0CompensationX/aflTimeConstant[3]" />
						  </value>
					  </userParameterDouble>
					</xsl:if>

					<!-- B0 compensation constant Y -->
					<xsl:if test="siemens/MEAS/sGRADSPEC/sB0CompensationY/aflTimeConstant[1]">
					  <userParameterDouble>
						  <name>sB0CompensationY_Tau_0</name>
						  <value>
							  <xsl:value-of select="siemens/MEAS/sGRADSPEC/sB0CompensationY/aflTimeConstant[1]" />
						  </value>
					  </userParameterDouble>
					</xsl:if>
					
					
					<xsl:if test="siemens/MEAS/sGRADSPEC/sB0CompensationY/aflTimeConstant[2]">
					  <userParameterDouble>
						  <name>sB0CompensationY_Tau_1</name>
						  <value>
							  <xsl:value-of select="siemens/MEAS/sGRADSPEC/sB0CompensationY/aflTimeConstant[2]" />
						  </value>
					  </userParameterDouble>
					</xsl:if>
					
					
					<xsl:if test="siemens/MEAS/sGRADSPEC/sB0CompensationY/aflTimeConstant[3]">
					  <userParameterDouble>
						  <name>sB0CompensationY_Tau_2</name>
						  <value>
							  <xsl:value-of select="siemens/MEAS/sGRADSPEC/sB0CompensationY/aflTimeConstant[3]" />
						  </value>
					  </userParameterDouble>
					</xsl:if>
					
					
					<!-- B0 compensation constant Z -->
					<xsl:if test="siemens/MEAS/sGRADSPEC/sB0CompensationZ/aflTimeConstant[1]">
					  <userParameterDouble>
						  <name>sB0CompensationZ_Tau_0</name>
						  <value>
							  <xsl:value-of select="siemens/MEAS/sGRADSPEC/sB0CompensationZ/aflTimeConstant[1]" />
						  </value>
					  </userParameterDouble>
					</xsl:if>
					
					
					<xsl:if test="siemens/MEAS/sGRADSPEC/sB0CompensationZ/aflTimeConstant[2]">
					  <userParameterDouble>
						  <name>sB0CompensationZ_Tau_1</name>
						  <value>
							  <xsl:value-of select="siemens/MEAS/sGRADSPEC/sB0CompensationZ/aflTimeConstant[2]" />
						  </value>
					  </userParameterDouble>
					</xsl:if>
					
					
					<xsl:if test="siemens/MEAS/sGRADSPEC/sB0CompensationZ/aflTimeConstant[3]">
					  <userParameterDouble>
						  <name>sB0CompensationZ_Tau_2</name>
						  <value>
							  <xsl:value-of select="siemens/MEAS/sGRADSPEC/sB0CompensationZ/aflTimeConstant[3]" />
						  </value>
					  </userParameterDouble>
					</xsl:if>
					<!-- ECC end -->
								
					
					
			</userParameters>
					

        </ismrmrdHeader>
    </xsl:template>

</xsl:stylesheet>
