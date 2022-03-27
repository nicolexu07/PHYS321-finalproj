import wget
import os
import requests
import pandas as pd

cwd = os.getcwd()



def get_links(site):
    """ (str) -> (list)
    Returns a list of strings, each entry in the list is a http:// link found on the website site.
    Ignores links beginning with http://irsa as these do not link to .tbl files.
    """
    site_str = str(requests.get(site).content) #contents of site as string
    links = []

    i = 0
    while i < len(site_str)-8:
        if site_str[i:i+7]=='http://' and site_str[i+7:i+7+4]!='irsa':
            link = 'http://'
            j = 0
            while site_str[i+7+j] != ' ':
                link += site_str[i+7+j]
                j += 1
            links.append(link)
            i += j
        i += 1
    return links



def download_from_link(link):
    """ (str) -> ()
    Downloads .tbl files found through the links at:
    https://exoplanetarchive.ipac.caltech.edu/bulk_data_download/wget_RADIAL.bat
    Saves them to a folder in the current working directory called data.
    """
    wget.download(link, cwd+'\data')



def download_radial():
    """ () -> ()
    Downloads all the .tbl files linked to in the exoplanetarchive.
    """
    links = get_links('https://exoplanetarchive.ipac.caltech.edu/bulk_data_download/wget_RADIAL.bat')
    for link in links:
        download_from_link(link)



def get_data(filename):
    """ (str) -> (pd.DataFrame)
    Takes a file name as input 
    Returns a pandas dataframe with 3 cols corresponding
    to days, radial velocity, and uncertainty
    """
    
    df = pd.read_table(f'./data/{filename}', header=None, skiprows=22, delim_whitespace=True)
    
    return df
        
    
 
 def get_obs_info(filename):
    """ (str) -> (np.array)
    Takes file name as input and returns numpy array with
    observation information as entries 
    """
    df = pd.read_table(f'./data/{filename}', header=None, nrows=19, delim_whitespace=True)
    
    # converting into numpy array with relevant info as entries
    temp = pd.DataFrame.to_numpy(df)
    
    info = []
    for entry in temp:
        info.append(entry[2])
    
    info = np.array(info)
    
    return info



if __name__ == "__main__":
    import doctest
    doctest.testmod()