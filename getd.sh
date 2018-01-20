#!/bin/bash

url_p="http://sha.ipac.caltech.edu/applications/Spitzer/SHA/servlet/SaveAsIpacTable?edu.caltech.ipac.firefly.data.Request=id%3DaorByProgramID%26RequestClass%3DServerRequest%26DoSearch%3Dtrue%26SearchByProgram.field.program%3D"
url_t="%26MoreOptions.field.prodtype%3Daor%26InstrumentPanel.field.irac%3D_none_%26InstrumentPanel.field.mips%3D_none_%26InstrumentPanel.field.irs%3Dstare%2Cmap%2Chi10%2Chi19%2Clow5%2Clow7%2Clow14%2Clow20%2Cblue%2Cred%26InstrumentPanel.field.panel%3Dinstrument%26inclCols%3Dtargetname%2Craj2000%2Cdecj2000%2Cnaifid%2Cmodedisplayname%2Creqkey%2Creqtitle%2Creqbegintime%2Creqendtime%2Cprogid%2Cpi%26pageSize%3D2147483647&file_name=aorByProgramID"
for i in $(cat query.txt); do
	echo $i
#	echo $url_p$i$url_t
	content="$(curl -s "$url_p$i$url_t")"
	echo "$content" >> output$i.txt
done
