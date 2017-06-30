"""
 Purpose:
  Download requested ISD files from the NOAA FTP site.
 Usage:
  python get_isd_from_noaa.py
 Author: PRI
 Date: June 2017
"""
from configobj import ConfigObj
import ftplib
import logging
import os

logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%H:%M:%S', level=logging.DEBUG)
# read the control file
cf_name = "../controlfiles/ISD/isd2.txt"
logging.info(" Getting config file contents")
cf = ConfigObj(cf_name)
site_list = list(cf["Sites"].keys())
isd_base_path = cf["Files"]["isd_base_path"]
# parse the control file to get a dictionary containing lists of the ISD sites
# required for each year
ftp_site_list = {}
for site in site_list:
    isd_site_list = cf["Sites"][site]["isd_sites"]
    start_year = int(cf["Sites"][site]["start_year"])
    end_year = int(cf["Sites"][site]["end_year"])
    year_list = [str(yr) for yr in range(start_year,end_year+1)]
    for year in year_list:
        if year not in ftp_site_list:
            ftp_site_list[year] = []
        ftp_site_list[year] += isd_site_list
# remove duplicate ISD site entries
year_list = sorted(ftp_site_list.keys())
for year in year_list:
    ftp_site_list[year] = list(set(ftp_site_list[year]))
# connect to the NOAA ISD ftp site
logging.info(" Opening NOAA FTP site")
ftp_base_path = os.path.join("pub","data","noaa")
ftp = ftplib.FTP("ftp.ncdc.noaa.gov","ftp","pisaac.ozflux@gmail.com")
# loop over the years and get the ISD files
for year in year_list:
    logging.info("  Processing year "+year)
    # get the directory for this year on the FTP server
    ftp_year_dir = os.path.join(ftp_base_path,year)
    # get a list of site files in this directory
    logging.info("   Getting the ISD file list for "+year)
    ftp_file_list = ftp.nlst(ftp_year_dir)
    # get the path where site files from this year will be written
    out_year_path = os.path.join(isd_base_path,year)
    # create the directory if it doesn't exist
    if not os.path.exists(out_year_path):
        os.makedirs(out_year_path)
    # now loop over the ISD sites required for this year
    for isd_site in ftp_site_list[year]:
        # construct the ISD file name
        ftp_file_name = isd_site+"-"+year+".gz"
        # and get the full path
        ftp_file_path = os.path.join(ftp_year_dir,ftp_file_name)
        # check that the site file exists in the FTP directory
        if ftp_file_path in ftp_file_list:
            # if it does, download it
            logging.info("   Getting data for "+isd_site+" for year "+year)
            out_file_path = os.path.join(out_year_path,ftp_file_name)
            logging.info("   Writing to "+out_file_path)
            ftp.retrbinary("RETR " + ftp_file_path, open(out_file_path, 'wb').write)
        else:
            logging.info("   ISD site "+isd_site+" not found for year "+year)
# say goodbye
ftp.quit()