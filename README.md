# Pubmed Search
 
This code contains functions for using the NIH's Entrez API to search the PubMed database and retrieve/download search results. 

## get_results

This function retrieves raw results from a PubMed search given an input string to search. Specify the search type using the `parameter_type` variable, which is a string. You can choose between `'keyword'`, `'identifier'`, `'author'`, or `'advanced'`.

For the first three options, use the `search_term` variable and enter the desired search term as a string. For the advanced search, enter the search query from the [PubMed Advanced Search Builder](https://pubmed.ncbi.nlm.nih.gov/advanced/) as a string. 

The `results_cap` variable dictates how many results will be returned (the default is 10). If the search returns a list of results longer than the value of `results_cap`, a warning will be printed.

The output of this function is a dictionary containing the raw search results from PubMed. 

## get_pmids

This function takes in the raw results returned by `get_results()` and returns the PMIDs as a list of strings. 

## export_pmids

This function takes in the raw results returned by `get_results()` and exports the list of PMIDs in .txt or .csv format. You can dictate the desired file type by specifying the `filetype` variable as `'.txt'` or `'.csv'`. You may also assign the value of the `filename` variable as a string to specify the name of the output file, as well as the desired saving location by adding a file path.

##export_results

This function takes in the raw results returned by `get_results()` and exports the full details in .txt or .csv format. You can dictate the desired file type by specifying the `filetype` variable as `'.txt'` or `'.csv'`. You may also assign the value of the `filename` variable as a string to specify the name of the output file, as well as the desired saving location by adding a file path.
