import wget
import os

cwd = os.getcwd()

def get_links(site):
    """ (str) -> (list)
    Returns a list of strings, each entry in the list is a http:// link found on the website site.
    """
    return None

def download_fron_link(link):
    """ (str) -> ()
    Downloads .tbl files found through the links at:
    https://exoplanetarchive.ipac.caltech.edu/bulk_data_download/wget_RADIAL.bat
    Saves them to a folder in the current working directory called data.
    """
    wget.download(link, cwd+'\data')