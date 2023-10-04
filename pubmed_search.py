from Bio import Entrez

def get_results(parameter_type, search_term, results_cap=10):
    '''Searches PubMed and returns raw search results indluding the list of PMIDs. 
    parameter_type: "keyword", "identifier", "author", "advanced" (copy search string from PubMed advanced search builder)
    Enter search term as string. Defauly results cap is 10 unless specified.
    '''
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
    
def get_pmids(results):
    '''Extracts the PMIDs from the raw results and returns as a list of strings.'''
    pmids = results["IdList"]
    return pmids

def export_pmids(results, filetype, filename = ''):
    '''Export the PMIDs as a list in .csv or .txt format. Enter the filetype as '.txt' or '.csv'. '''
    pmids = results['IdList']
    with open('pmids' + filename + filetype, 'w') as f:
        for line in pmids:
            f.write(line)
            f.write('\n')
    
def export_results(results, filetype, filename = ''):
    '''Exports the full search results in either .csv or .txt format. Enter the filetype as '.txt' or '.csv.'''
    
    with open('search_results' + filename + filetype, 'w') as f:
        f.write('Count \t \t \t '+ results['Count'] + '\n')
        f.write('RetMax \t \t \t '+ results['RetMax'] + '\n')
        f.write('RetStart \t \t '+ results['RetStart'] + '\n')
        f.write('TranslationSet \t \t'+ str(results['TranslationSet']) + '\n')
        f.write('QueryTranslation \t '+ str(results['QueryTranslation']) + '\n')
        f.write('PMIDs \n \n')
        for pmid in results['IdList']:
            f.write(pmid + '\n')


