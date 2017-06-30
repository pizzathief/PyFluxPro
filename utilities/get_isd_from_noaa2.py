from configobj import ConfigObj
import ftplib
import logging
import os

logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%H:%M:%S', level=logging.DEBUG)

cf_name = "../controlfiles/ISD/isd2.txt"
logging.info(" Getting config file contents")
cf = ConfigObj(cf_name)
site_list = list(cf["Sites"].keys())
isd_base_path = cf["Files"]["isd_base_path"]

logging.info(" Opening NOAA FTP site")
ftp_base_path = os.path.join("pub","data","noaa")
ftp = ftplib.FTP("ftp.ncdc.noaa.gov","ftp","pisaac.ozflux@gmail.com")

got_isd = {}
ftp_file_list = {}
for site in site_list:
    logging.info(" Processing site "+site)
    start_year = int(cf["Sites"][site]["start_year"])
    end_year = int(cf["Sites"][site]["end_year"])
    year_list = [str(year) for year in range(start_year, end_year+1)]
    isd_site_list = cf["Sites"][site]["isd_sites"]
    for year in year_list:
        logging.info("  Processing year "+year)
        got_isd[year] = []
        ftp_year_dir = os.path.join(ftp_base_path,year)
        if year not in list(ftp_file_list.keys()):
            logging.info("   Getting the ISD file list for "+year)
            ftp_file_list[year] = ftp.nlst(ftp_year_dir)
        out_year_path = os.path.join(isd_base_path,year)
        if not os.path.exists(out_year_path):
            os.makedirs(out_year_path)
        for isd_site in isd_site_list:
            ftp_file_name = isd_site+"-"+year+".gz"
            ftp_file_path = os.path.join(ftp_year_dir,ftp_file_name)
            if ftp_file_path in ftp_file_list[year]:
                if isd_site not in got_isd[year]:
                    logging.info("   Getting data for "+isd_site+" for year "+year)
                    out_file_path = os.path.join(out_year_path,ftp_file_name)
                    logging.info("   Writing to "+out_file_path)
                    ftp.retrbinary("RETR " + ftp_file_path, open(out_file_path, 'wb').write)
                    got_isd[year].append(isd_site)
                else:
                    logging.info("   Already got ISD site "+isd_site+" for year "+year)
            else:
                logging.info("   ISD site "+isd_site+" not found for year "+year)

ftp.quit()