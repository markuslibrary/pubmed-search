from Bio import Entrez
from Bio import Medline

# Remember to define your Entrez email like so:
# Entrez.email = 'your@email.com'


# Function to search PubMed and return a raw list of results
def get_results(parameter_type, search_term, results_cap=10):
    '''
    parameter_type: 'keyword', 'identifier', 'author', 'advanced' (copy search string from PubMed advanced search builder)
    Enter search term as string. Default results cap is 10 unless specified.
    '''
    
    # Dictionary with allowed search parameters
    parameter_dict = {'keyword' : search_term, 
                      'identifier' : search_term + '[SI]', 
                      'author' : search_term + '[Author]',
                      'advanced': search_term}
    
    search = parameter_dict[parameter_type]
    
    try:
        # Create search query and get results
        handle = Entrez.esearch(db = 'pubmed', term=search, retmax = results_cap)
        results = Entrez.read(handle)
        
        # Print warning if the cap on results is smaller than the number of results returned
        total_results = results['Count']
        if int(total_results) > results_cap:
            print('Warning: the complete list of results is not being returned.')
        
        return results
        
    except Exception as e:
        print(f'An error occurred: {str(e)}.')

# Function to extract the PMIDs from the raw results    
def get_pmids(results):
    pmids = results["IdList"]
    return pmids

# Function to export PMIDs from raw search results in .csv or .txt format
def export_pmids(results, filetype, filename = 'pmids'):
    '''Enter the filetype as '.txt' or '.csv'. Enter the desired output file name, otherwise the defauls will be used'''
    pmids = results['IdList']
    with open('pmids' + filename + filetype, 'w') as f:
        for line in pmids:
            f.write(line)
            f.write('\n')

# Function to export the full search results in either .csv or .txt format
def export_results(results, filetype, filename = 'results'):
    ''' Enter the filetype as '.txt' or '.csv.'''
    with open('search_results' + filename + filetype, 'w') as f:
        f.write('Count \t \t \t '+ results['Count'] + '\n')
        f.write('RetMax \t \t \t '+ results['RetMax'] + '\n')
        f.write('RetStart \t \t '+ results['RetStart'] + '\n')
        f.write('TranslationSet \t \t'+ str(results['TranslationSet']) + '\n')
        f.write('QueryTranslation \t '+ str(results['QueryTranslation']) + '\n')
        f.write('PMIDs \n \n')
        for pmid in results['IdList']:
            f.write(pmid + '\n')

# Function that takes in a list of PMIDs and writes a RIS citation file which can be opened with any citation management software
def export_citations(pmids, filename = 'citations'):
    ''' pmids should be a list of PMIDs, listed as strings. '''
    try:
        records = []
        for pmid in pmids:
            handle = Entrez.efetch(db = 'pubmed', id = pmid, rettype = 'medline', retmode = 'text')
            record = Medline.read(handle)
            records.append(record)

        with open(filename + '.ris', 'w', encoding='utf-8') as ris_file:
            for record in records:
                ris_file.write('TY  - JOUR\n')  # Type of reference (Journal article)
                ris_file.write(f"TI  - {record.get('TI', '')}\n")  # Title
                ris_file.write(f"AB  - {record.get('AB', '')}\n")  # Abstract
                ris_file.write(f"AU  - {', '.join(record.get('AU', []))}\n")  # Authors
                ris_file.write(f"PY  - {record.get('DP', '')}\n")  # Publication Year
                ris_file.write(f"PMID  - {record.get('PMID', '')}\n")  # PMID
                ris_file.write('ER  - \n\n')  # End of reference

    except Exception as e:
        print('An error occurred:', str(e))


