import logging
import os

def init_logger(logger_name="pfp_log"):
    """
    Purpose:
     Returns a logger object.
    Usage:
     logger = qclog.init_logger()
    Author: PRI with acknowledgement to James Cleverly
    Date: September 2016
    """
    logger = logging.getLogger(name=logger_name)
    logger.setLevel(logging.DEBUG)
    # create file handler
    #max_bytes = 1024 * 1024 * 2
    #fh = logging.handlers.RotatingFileHandler(os.path.join("logfiles", 'pfp.log'), mode="a", maxBytes=max_bytes, backupCount=1)
    fh = logging.FileHandler(os.path.join("logfiles", 'pfp.log'))
    fh.setLevel(logging.DEBUG)
    # create console handler
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    # create formatter and add to handlers
    formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s','%H:%M:%S')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    # add the handlers to the logger
    logger.addHandler(fh)
    logger.addHandler(ch)

    return logger
