get_solarmonitor_links_dirty
============================
This code "queries" `solarmonitor.org <https://www.solarmonitor.org>`_ for four different full sun images 
(Line-of-Sight (LoS) Magnetic Fields, Hinode-XRT, H :math:r`\alpha`, and ~171&Aring;).
The query works by making text requests for a specified date to solarmonitor then parsing the text for png files related to observations in all four wavelengths.
If the program finds a file, then it downloads the file locally.
In preceeding iterations, the program checks that the file does not already exist locally before trying to download the four images.
Normally, get_solarmonitor_links_dirty runs from within create_json_for_web, but may also run independently.
An example of the code is below.

sswidl>;Downloads all files not currently stored locally in output_dir/ and stores the path to all local files in links, even if they were previously downloaded.
sswidl>get_solarmonitor_links_dirty,'2019/01/04','output_dir/',links