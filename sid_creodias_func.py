import wget
import requests
from keycloak import KeycloakOpenID
import keycloak
import urllib
import os
import glob
from osgeo import ogr#, osr
import numpy as np
import pathlib

URL_DOWNLOAD = 'https://zipper.creodias.eu/download/UID?token='
URL_TOKEN = 'https://identity.cloudferro.com/auth/realms/Creodias-new/protocol/openid-connect/token'


def download_from_polygon(polygon, start_date, end_date, out_folder, aut2, start_hour = '00', end_hour = '23', coverage=.8):
    '''
    This function downloads Sentinel-1 GRD products from Creodias
    
    Inputs:
    -------
    polygon: array_like
        Lats and longs of the area of interest
        format [lat1,lat2,lon1,lon2]
    start_date: str_like
        format 'YYYY-MM-DD'
        The start date the observations of interest were made
    end_date: string
        format 'YYYY-MM-DD'
        The end date the observations of interest were made
    out_folder: string
        Path of the folder where to download the S1 products
    start_hour: string (optional)
        format 'HH'
        The start hour in the start_date the observations of interest were made
    end_hour: string (optional)
        format 'HH'
        The end hour in the end_date the observations of interest were made
        
        
    Outputs:
    -------
    r: bool
        0 if the donwload did not succeed due to errors of connection / no good overlap
        1 if the product is downloaded
    '''
    
    if start_hour != '00':
        end_hour = str(int(start_hour)+1)
        if len(end_hour) == 1:
            end_hour = f'0{end_hour}'
            
    path = os.path.dirname(__file__)
    #token = create_token()
    token = _get_token(aut2)    #get token with this new key from app

    # Remove temporary files
    for temp_files in glob.glob(f'{path}/*.tmp', recursive = True):
        os.remove(temp_files)

    # Extract the downlod url
    URL = f'https://finder.creodias.eu/resto/api/collections/Sentinel1/search.json?maxRecords=10&startDate={start_date}T{start_hour}%3A00%3A00Z&completionDate={end_date}T{end_hour}%3A59%3A59Z&productType=GRD&sensorMode=EW&geometry=POLYGON(({polygon[2]}+{polygon[1]}%2C{polygon[3]}+{polygon[1]}%2C{polygon[3]}+{polygon[0]}%2C{polygon[2]}+{polygon[0]}%2C{polygon[2]}+{polygon[1]}))&sortParam=startDate&sortOrder=descending&status=all&dataset=ESA-DATASET'

    r = find_products(URL)
        
    if "ErrorCode" in r.json():
        return 0
    
    if r.json()["properties"]["totalResults"] == 0:
        print("No product fits with your quiery")
        return 0

    # Calculate intersection of the found products and the polygon of interest    
    intersection=[]
    for i in range(r.json()["properties"]["totalResults"]):
        #print(i)
        try:
            inter = get_intersection(
                r.json()["features"][i]["geometry"]["coordinates"][0], polygon)
        except:
            #some scenes dont have coodinates...
            inter=0
        intersection.append(inter)
        
    print('Intersection fractions of individual scenes with your polygon: ',intersection)

    id = np.argmax(intersection)    #we only grab the best match!
    
    # If all products have less than 80% (default) of overlap abort 
    if intersection[id]<=coverage:
        return 0
    
    fname=r.json()["features"][id]["properties"]["title"]
    
    # Check if the output folder exist 
    if folder_exist(out_folder):
        # If the product exist abort
        if product_exist(out_folder, pathlib.Path(r.json()["features"][id]["properties"]["productIdentifier"]).stem):
            print('File already downloaded')
            return 1,fname
        
    # Extract the url for download
    #url_download = f'{r.json()["features"][id]["properties"]["services"]["download"]["url"]}?token={token["access_token"]}'
    url_download = f'{r.json()["features"][id]["properties"]["services"]["download"]["url"]}?token={token}'
    
    #The product is not downloaded from finder, but from zipper
    url_download = url_download.replace("finder", 'zipper')
    
    # Download file
    pathlib.Path(out_folder).mkdir(parents=True, exist_ok=True)
    r = download_products(url_download, out_folder)
    return 0 if r is None else 1,fname

def OpenID():
    try:
        r = KeycloakOpenID(server_url="https://auth.creodias.eu/auth/",
                    client_id="CLOUDFERRO_PUBLIC",
                    realm_name="DIAS")
    except keycloak.exceptions.KeycloakError:
        # Maybe set up for a retry, or continue in a retry loop
        print('Timeout')
        OpenID()
    except keycloak.exceptions.KeycloakGetError:
        # Maybe set up for a retry, or continue in a retry loop
        print('Timeout')
        OpenID()
    return r

def find_products(url: str):
    try:
        r = requests.get(url)
    except requests.exceptions.Timeout:
        # Maybe set up for a retry, or continue in a retry loop
        print('Timeout')
        r = find_products(url)
    except requests.exceptions.TooManyRedirects:
        return None
    except requests.exceptions.RequestException as e:
        # catastrophic error. bail.
        raise SystemExit(e) from e
    return r

def download_products(url, out):
    try:
        r = wget.download(url, out)
    except urllib.error.ContentTooShortError:
        # Maybe set up for a retry, or continue in a retry loop
        print('Content too short')
        r = download_products(url, out)
    except urllib.error.HTTPError:
        r = None
    return r

#similar to get_token - this is redundant now
def create_token():

    keycloak_openid = OpenID()
    print("Id opened")
    return keycloak_openid.token("polona.itkin@uit.no", "GGP2D3j2uDZczxK!")


def get_intersection(get_poly_im, polygon):
    if len(get_poly_im) <5:
        get_poly_im = get_poly_im[0]
    if len(get_poly_im) >5:
        poly_im = get_lat_long(np.asarray(get_poly_im))
        wkt_im = f"POLYGON (({str(poly_im[2])} {str(poly_im[1])}, {str(poly_im[3])} {str(poly_im[1])}, {str(poly_im[3])} {str(poly_im[0])}, {str(poly_im[2])} {str(poly_im[0])}, {str(poly_im[2])} {str(poly_im[1])}))"
    else:
        wkt_im = f"POLYGON (({str(get_poly_im[0][0])} {str(get_poly_im[0][1])}, {str(get_poly_im[1][0])} {str(get_poly_im[1][1])}, {str(get_poly_im[2][0])} {str(get_poly_im[2][1])}, {str(get_poly_im[3][0])} {str(get_poly_im[3][1])}, {str(get_poly_im[4][0])} {str(get_poly_im[4][1])}))"

    wkt_oil = f"POLYGON (({str(polygon[2])} {str(polygon[1])}, {str(polygon[3])} {str(polygon[1])}, {str(polygon[3])} {str(polygon[0])}, {str(polygon[2])} {str(polygon[0])}, {str(polygon[2])} {str(polygon[1])}))"
    
    poly_oil = ogr.CreateGeometryFromWkt(wkt_oil)
    poly_im = ogr.CreateGeometryFromWkt(wkt_im)

    intersection = poly_oil.Intersection(poly_im)
    if intersection is None: 
        return 0
    else: return intersection.GetArea()/poly_oil.GetArea()

def folder_exist(folder_name):
    return bool(os.path.exists(f'{folder_name}'))

def product_exist(folder_name, product_name):
    return bool(os.path.exists(f'{folder_name}/{product_name}.SAFE') or os.path.exists(f'{folder_name}/{product_name}.zip'))

def get_lat_long(array):
    if len(array.shape) == 2:
        return np.min(array[:,1]), np.max(array[:,1]), np.min(array[:,0]), np.max(array[:,0])
    else: return np.min(array[0,0,:,1]), np.max(array[0,0,:,1]), np.min(array[0,0,:,0]), np.max(array[0,0,:,0])   
    
#the functions below are copied from the https://source.coderefinery.org/cirfa/sat_search_download/
#and are written by Henrik Fisser  - adopted to work with arguments to script
def _get_token(aut2, username="polona.itkin@uit.no", password="GGP2D3j2uDZczxK!"):
    refresh_token = _get_refresh_token(username)  # get refresh token from previous call stored in file
    token_data = {
        'client_id': 'CLOUDFERRO_PUBLIC',
        'username': username,
        'password': password,
        'grant_type': 'refresh_token',
        'refresh_token': refresh_token
    }
    if len(refresh_token) > 0:  # don't bother API if we know it isn't working because no refresh token
        response = requests.post(URL_TOKEN, data=token_data).json()
    else:
        response = dict()
    try:
        token = response['access_token']
    except KeyError:  # we didn't find a refresh token or it's expired
        token_data['grant_type'] = 'password'  # renew login with password
        token_data['totp'] = input(f'****** Provide the 2-factor authentication token for user "{username}":   ')  # ask user for token from App on phone
        #token_data['totp'] = aut2
        del token_data['refresh_token']  #  the disfunctional refresh token should not be provided
        response = requests.post(URL_TOKEN, data=token_data).json()  # get token
        try:
            token = response['access_token']  # it should be there
        except KeyError:
            raise RuntimeError(f'Unable to get token. Response from CREODIAS was {response}')  # could be wrong user credentials, or service is down etc.
    finally:
        _update_refresh_token(response['refresh_token'], username)  # response contains a refresh token, we keep it to maximize the time we don't ask the user for a 2-factor authentication token
    #print(token)
    #exit()
    return token


def _update_refresh_token(refresh_token, username):
    os.environ[f'REFRESH_TOKEN_{username}'] = refresh_token


def _get_refresh_token(username):
    try:
        return os.environ[f'REFRESH_TOKEN_{username}']
    except KeyError:
        return ''  

if __name__ == '__main__':
    dir_out = os.path.dirname(__file__)
    username, password = os.environ["CREO_USER"], os.environ["CREO_PASSWORD"]
    file_zip = os.path.join(dir_out, "test_to_be_deleted.zip")
    kwargs = dict(username=username, password=password)
    token = _get_token(**kwargs)  # test get token
    token = _get_token(**kwargs)  # retry test to see if it works with refresh token
    download(uid="a8eebdc2-31bb-57fa-8f96-3aeae4f2cdaf", username=username, password=password, outfile=file_zip)  # test if download also works
    os.remove(file_zip)  # delete downloaded data
