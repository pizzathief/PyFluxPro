import ftplib
import os

isd_sites = {}
isd_sites["2007"] = ["999999-03048","999999-03033","723650-23050","722677-03027","723620-99999",
                     "999999-03062","723654-93091","723656-23049"]
isd_sites["2008"] = ["999999-03048","999999-03033","723650-23050","722677-03027","723620-99999",
                     "999999-03062","723654-93091","723656-23049"]
isd_sites["2009"] = ["999999-03048","999999-03033","723650-23050","722677-03027","723620-99999",
                     "999999-03062","723654-93091","723656-23049"]
isd_sites["2010"] = ["999999-03048","999999-03033","723650-23050","722677-03027","723620-99999",
                     "999999-03062","723654-93091","723656-23049"]
isd_sites["2011"] = ["999999-03048","999999-03033","723650-23050","722677-03027","723620-99999",
                     "999999-03062","723654-93091","723656-23049"]
isd_sites["2012"] = ["999999-03048","999999-03033","723650-23050","722677-03027","723620-99999",
                     "999999-03062","723654-93091","723656-23049"]
isd_sites["2013"] = ["999999-03048","999999-03033","723650-23050","722677-03027","723620-99999",
                     "999999-03062","723654-93091","723656-23049"]
isd_sites["2014"] = ["999999-03048","999999-03033","723650-23050","722677-03027","723620-99999",
                     "999999-03062","723654-93091","723656-23049"]
isd_sites["2015"] = ["999999-03048","999999-03033","723650-23050","722677-03027","723620-99999",
                     "999999-03062","723654-93091","723656-23049"]
isd_sites["2016"] = ["999999-03048","999999-03033","723650-23050","722677-03027","723620-99999",
                     "999999-03062","723654-93091","723656-23049"]

year_list = list(isd_sites.keys())
isd_base_path = os.path.join("pub","data","noaa")
out_base_path = "/home/peter/AmeriFlux/ISD"

ftp = ftplib.FTP("ftp.ncdc.noaa.gov","ftp","pisaac.ozflux@gmail.com")

for year in year_list:
    isd_year_dir = os.path.join(isd_base_path,year)
    isd_file_list = ftp.nlst(isd_year_dir)
    for isd_site in isd_sites[year]:
        isd_file_name = isd_site+"-"+year+".gz"
        isd_file_path = os.path.join(isd_year_dir,isd_file_name)
        #print isd_file_path
        if isd_file_path in isd_file_list:
            print "Getting data for "+isd_site+" for year "+year
            out_file_path = os.path.join(out_base_path,year,isd_file_name)
            #print "Asking for: ",isd_file_name
            #print "Writing to: ",out_file_path
            ftp.retrbinary("RETR " + isd_file_path, open(out_file_path, 'wb').write)

ftp.quit()