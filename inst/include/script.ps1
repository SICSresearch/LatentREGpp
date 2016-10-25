$Files = Get-Content dlibUsedHeaders.txt

foreach ($File in $Files) {
	$Temp = "output\" + (Split-Path -Path $File)
	echo $File
	echo $Temp
	New-Item $Temp -type directory -Force
	Copy-Item $File $Temp
}