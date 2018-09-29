## jgi_metagenome_scraping

from selenium import webdriver
import time
import os
import re
import json
from bs4 import BeautifulSoup        

def activate_driver():
    """
    Activate chrome driver used to automate webpage navigation (see: https://sites.google.com/a/chromium.org/chromedriver/)
    The chrome driver .exec file must be in the home directory

    returns: driver [object]
    """
    homedir = os.path.expanduser('~')
    return webdriver.Chrome(homedir+'/chromedriver')

def get_ecosystemclass_url_from_jgi_img_homepage(driver,homepage_url,ecosystemClass,database='jgi'):
    """
    load homepage_url -> retrieve ecosytemClass_url

    driver: the chrome driver object
    homepage_url: url of the jgi homepage. should be 'https://img.jgi.doe.gov/cgi-bin/m/main.cgi' as of 6/15/2017
    ecosystemClass: can be 'Engineered', 'Environmental', or 'Host-associated' (these are 3 different links on the homepage)
    database: choose to use only the jgi database, or all database [default=jgi]
    """

    driver.get(homepage_url)
    time.sleep(5)
    htmlSource = driver.page_source

    ## All ampersands (&) must be followed by 'amp;'
    regex = r'href=\"main\.cgi(\?section=TaxonList&amp;domain=Metagenome&amp;seq_center=%s&amp;page=metaCatList&amp;phylum=%s)\"'%(database,ecosystemClass)

    match = re.search(regex, htmlSource)
    ecosystemClass_suffix = match.group(1)
    ecosystemClass_url = homepage_url+ecosystemClass_suffix
    # ecosystemClass_urls.append(ecosystemClass_url)

    return ecosystemClass_url

def get_ecosystemclass_json_from_ecosystem_class_url(driver,ecosystemClass_url):
    """
    load ecosystemClass_url- > retrieve ecosystemClass_json_url -> load ecosystemClass_json_url -> retrieve ecosystemClass_json
    
    driver: the chrome driver object
    ecosystemClass_url: url of either the 'Engineered', 'Environmental', or 'Host-associated' ecosystemClasses (these are 3 different links on the homepage)
    """

    driver.get(ecosystemClass_url)
    time.sleep(5)
    htmlSource = driver.page_source
    # driver.quit()

    regex = r'var myDataSource = new YAHOO\.util\.DataSource\(\"(.*)\"\);'
    match = re.search(regex, htmlSource)
    ecosystemClass_json_suffix = match.group(1)
    ecosystemClass_url_prefix = ecosystemClass_url.split('main.cgi')[0]
    ecosystemClass_json_url = ecosystemClass_url_prefix+ecosystemClass_json_suffix

    driver.get(ecosystemClass_json_url)
    time.sleep(5)
    # jsonSource = driver.page_source

    ## convert the jsonSource into a dict of dicts here
    ecosystemClass_json = json.loads(driver.find_element_by_tag_name('body').text)

    return ecosystemClass_json


def get_metagenome_urls_from_ecosystemclass_json(driver,homepage_url,ecosystemClass_json):
    """
    parse ecosystemClass_json -> retrieve metagenome_urls

    driver: the chrome driver object
    homepage_url: url of the jgi homepage. should be 'https://img.jgi.doe.gov/cgi-bin/m/main.cgi' as of 6/15/2017
    ecosystemClass_json: json text of either the 'Engineered', 'Environmental', or 'Host-associated' ecosystemClasses
    """

    all_GenomeNameSampleNameDisp =  [d['GenomeNameSampleNameDisp'] for d in ecosystemClass_json['records']]

    metagenome_urls = list()

    for htmlandjunk in all_GenomeNameSampleNameDisp:
        regex = r"<a href='main\.cgi(.*)'>"
        match = re.search(regex, htmlandjunk)
        html_suffix = match.group(1)
        full_url = homepage_url+html_suffix
        metagenome_urls.append(full_url)

    return metagenome_urls

def get_metagenome_htmlSource_and_metadata(driver, metagenome_url):
    """
    load metagenome_url -> retrieve metagenome_htmlSource & metadata

    driver: the chrome driver object
    metagenome_url: url for an single metagenome
    """

    driver.get(metagenome_url)
    time.sleep(5)
    metagenome_htmlSource = driver.page_source

    metadata_table_dict = get_metagenome_metadata_while_on_metagenome_page(metagenome_htmlSource)

    return metagenome_htmlSource, metadata_table_dict


def get_enzyme_url_from_metagenome_url(metagenome_url, metagenome_htmlSource, datatype):

    """
    metagenome_htmlSource -> parse out enzyme_url

    driver: the chrome driver object
    metagenome_url: url for an single metagenome
    data_type: can be assembled, unassembled, or both (refers to whether the metagenomes are assembled or not)
    """

    
    regex = r'<a href=\"(main\.cgi\?section=MetaDetail&amp;page=enzymes.*data_type=%s.*)\" onclick'%datatype
    match = re.search(regex, metagenome_htmlSource)

    if match:

        print "Metagenome url: %s ...\n...has datatype: %s"%(metagenome_url,datatype)

        enzyme_url_suffix = match.group(1)
        enzyme_url_prefix = metagenome_url.split('main.cgi')[0]
        enzyme_url = enzyme_url_prefix+enzyme_url_suffix

    else:

        print "Metagenome url: %s ...\n...does not have datatype: %s"%(metagenome_url,datatype) 

        enzyme_url = None

    

    return enzyme_url


def get_metagenome_metadata_while_on_metagenome_page(htmlSource):
    """
    htmlSource -> dictionary of metagenome metadata

    htmlSource: the metagenome_url driver's .page_source
    """

    # return dict of metagenome table data
    bs = BeautifulSoup(htmlSource,"html.parser")
    metadata_table = bs.findAll('table')[0]

    metadata_table_dict = dict()
    for row in metadata_table.findAll('tr'):

        if (len(row.findAll('th')) == 1) and (len(row.findAll('td')) == 1):

            row_key = row.findAll('th')[0].text.rstrip()
            row_value = row.findAll('td')[0].text.rstrip() if row.findAll('td')[0] else None
            metadata_table_dict[row_key] = row_value

    metadata_table_dict.pop('Geographical Map', None)

    ## metadata_table_dict['Taxon Object ID'] should be the way we identify a metagenome

    return metadata_table_dict

def get_enzyme_json_from_enzyme_url(driver,enzyme_url):
    """
    load enzyme_url -> retrieve enzyme_json_url -> load enzyme_json_url -> retrieve enzyme_json

    driver: the chrome driver object
    enzyme_url: url for an single enzyme type from an single metagenome
    """

    driver.get(enzyme_url)
    time.sleep(5)
    htmlSource = driver.page_source
    # driver.quit()

    regex = r'var myDataSource = new YAHOO\.util\.DataSource\(\"(.*)\"\);'
    match = re.search(regex, htmlSource)
    enzyme_json_suffix = match.group(1)
    enzyme_url_prefix = enzyme_url.split('main.cgi')[0]
    enzyme_json_url = enzyme_url_prefix+enzyme_json_suffix

    driver.get(enzyme_json_url)
    time.sleep(5)
    # jsonSource = driver.page_source

    ## convert the jsonSource into a dict of dicts here
    enzyme_json = json.loads(driver.find_element_by_tag_name('body').text)

    return enzyme_json

def parse_enzyme_info_from_enzyme_json(enzyme_json):

    enzyme_dict = dict() # Dictionary of ec:[enzymeName,genecount] for all ecs in a single metagenome

    for i, singleEnzymeDict in enumerate(enzyme_json['records']):
        ec = singleEnzymeDict['EnzymeID']
        enzymeName = singleEnzymeDict['EnzymeName']
        genecount = singleEnzymeDict['GeneCount']

        enzyme_dict[ec] = [enzymeName,genecount]

    return enzyme_dict



def main():

    driver = activate_driver()
    homepage_url = 'https://img.jgi.doe.gov/cgi-bin/m/main.cgi'
    database = 'jgi'
    ecosystemClasses = ['Engineered', 'Environmental', 'Host-associated']
    datatypes = ['assembled','unassembled','both']
    save_dir = 'jgi_metagenome_jsons'

    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    jgi_metagenomes = list()
    for ecosystemClass in ecosystemClasses:

        print "Scraping all metagenomes from ecosystemClass: %s ..."%ecosystemClass

        ecosystemClass_url = get_ecosystemclass_url_from_jgi_img_homepage(driver,homepage_url,ecosystemClass,database=database)

        ecosystemClass_json = get_ecosystemclass_json_from_ecosystem_class_url(driver,ecosystemClass_url)

        metagenome_urls = get_metagenome_urls_from_ecosystemclass_json(driver,homepage_url,ecosystemClass_json)

        for metagenome_url in metagenome_urls:

            print "Scraping metagenome: %s ..."%metagenome_url

            metagenome_htmlSource, metadata_table_dict = get_metagenome_htmlSource_and_metadata(driver, metagenome_url)

            single_metagenome_dict = {'metadata':metadata_table_dict}

            taxon_object_id = metadata_table_dict['Taxon Object ID']
            
            for datatype in datatypes:

                enzyme_url = get_enzyme_url_from_metagenome_url(metagenome_url, metagenome_htmlSource, datatype)

                if enzyme_url:

                    enzyme_json = get_enzyme_json_from_enzyme_url(driver,enzyme_url)

                    enzyme_dict = parse_enzyme_info_from_enzyme_json(enzyme_json)

                    single_metagenome_dict[datatype] = enzyme_dict

            jgi_metagenomes.append(single_metagenome_dict)

            with open(save_dir+'/'+taxon_object_id+'.json', 'w') as outfile:
        
                json.dump(single_metagenome_dict,outfile)

            print "Done scraping metagenome."
            print "-"*80

        print "Done scraping metagenomes from ecosystemClass: %s."%ecosystemClass
        print "="*90



    print "Writing json to file..."

    large_json_id = "_".join(ecosystemClasses)
    with open('jgi_metagenomes_'+large_json_id+'.json', 'w') as outfile:
        
        json.dump(jgi_metagenomes,outfile)

    print "Done."

    ## Can i write it so that it scrapes many at a time?

if __name__ == '__main__':
    main()




