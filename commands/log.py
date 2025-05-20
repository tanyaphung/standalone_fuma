import logging
import os


def log(filedir): 
    # create and configure logger
    logging.basicConfig(filename=os.path.join(filedir, "params_logs.txt"),
                        format='%(asctime)s %(levelname)s %(message)s',
                        filemode='w')
    
    logger = logging.getLogger() #creating an object
    logger.setLevel(logging.DEBUG) #setting the threshold of logger to DEBUG
    return logger