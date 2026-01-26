# study metadata sheet https://github.com/LieberInstitute/spatial_NAc/blob/main/raw-data/sample_key_spatial_NAc.csv
# GEO metadata sheet -> use that https://ftp.ncbi.nlm.nih.gov/geo/series/GSE307nnn/GSE307586/miniml/GSE307586_family.xml.tgz
import requests
import os
import tarfile
import pandas as pd

fname = "data/GSE307586_family.xml"

if not os.path.exists(fname):
    url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE307nnn/GSE307586/miniml/GSE307586_family.xml.tgz"
    r = requests.get(url)
    with open("data/GSE307586_family.xml.tgz", "wb") as f:
        f.write(r.content)
    with tarfile.open("data/GSE307586_family.xml.tgz", "r:gz") as tar:
        tar.extractall(path="data/")

import xml.etree.ElementTree as ET
namespace = {'ns': 'http://www.ncbi.nlm.nih.gov/geo/info/MINiML'}
tree = ET.parse(fname)
root = tree.getroot()

samples = root.findall('.//ns:Sample', namespace)

#for each sample, parse the attributes into a dict
#iid, Title, and Characteristics with tags: tissue, subject status, sample group
sdict = {}
for s in samples:
    iid = s.find('ns:Accession', namespace).text
    title = s.find('ns:Title', namespace).text
    tissue = s.find(".//ns:Characteristics[@tag='tissue']", namespace).text
    subject_status = s.find(".//ns:Characteristics[@tag='subject status']", namespace).text
    sdict[iid.strip()] = {
        'brain': title.split("_")[1].strip().replace("-Left",''),
        'tissue': tissue.strip(),
        'region':title.split('-')[-1].strip(),
        'subject_status': subject_status.strip()
    }

samplesheet = pd.DataFrame.from_dict(sdict, orient='index')
samplesheet.insert(0, 'data_directory', 'data/GSE307586/samples/'+samplesheet.index)
samplesheet.insert(0, 'sample', samplesheet.index)

# get additional metadata from study github
metadata_url = "https://raw.githubusercontent.com/LieberInstitute/spatial_NAc/main/raw-data/sample_key_spatial_NAc.csv"
metadata = pd.read_csv(metadata_url)
metadata.columns = [col.strip().lower() for col in metadata.columns]
metadata['region'] = [x.split("_")[-1] for x in metadata['slide']]
# drop duplicate columns in metadata
selected_metadata = metadata.groupby('brain').agg({'age':'mean','sex':'first'})

# merge metadata into samplesheet
samplesheet = samplesheet.merge(selected_metadata, on=['brain'], how='left')

# save to root for running STAPLE
samplesheet.to_csv("samplesheet.csv", index=False)