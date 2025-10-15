# R Shiny App
Development of a desktop application for 'Gene essentiality visualisation at the domain level'
App allows insertion of Tn-seq data and Pfam/InterPro domain annotations
Data needs to be preprocessed to enter clean data into the app, alinging gene IDs and excluding unneccessary information
Results include a table with per gene totals x per domain rows
Results include a plot with domain overlays
Results can be exported as an essentiality report CSV file or an essential domains BED file

'tn-seq.csv' file obtained from Mycobacterium tuberculosis H37Rv studies (DeJesus et al., 2017), the file was formatted as a wiggle file in TRANSIT.
This file is then put in the 'Processing data' code to get 'tnseq_clean.csv'.

'uniprot-proteome_UP000001584.tsv' file obtained from Uniprot Mycobacterium tuberculosis H37Rv, where the file is put in the 'Processing data' code to get 'domains_final.csv'.
