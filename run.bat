@REM scaffold hopper launcher for windows

@echo off
if "%HOPPER_HOME%"=="" set "HOPPER_HOME=%~dp0\\."


set "LIB=lib"

set "APP_CP=%LIB%\colt.jar;%LIB%\jackson-core-asl-1.9.2.jar;%LIB%\jackson-mapper-asl-1.9.2.jar;%LIB%\jchem.jar;%LIB%\jnlp.jar;%LIB%\jcommon-1.0.17.jar;%LIB%\jfreechart-1.0.14.jar;%LIB%\jgoodies-common-1.7.0.jar;%LIB%\jgoodies-looks-2.5.3.jar;%LIB%\jide-components.jar;%LIB%\jide-oss-3.5.9.jar;%LIB%\jxlayer.jar;%LIB%\prefuse.jar;%LIB%\swingworker-0.8.0.jar;%LIB%\swingx.jar"



@REM assume java is already installed

@echo "%APP_CP%"
java -cp "%APP_CP%;rgroup.jar" gov.nih.ncgc.rgroup.RGroupTool %*
