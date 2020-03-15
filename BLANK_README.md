
<!-- PROJECT LOGO -->
<br />
<p align="center">
  <a href="https://github.com/github_username/repo">
    <img src="images/logo.png" alt="Logo" width="80" height="80">
  </a>

  <h3 align="center">BE-FF</h3>

  <p align="center">
    Web-based tool that receives SNV data and matches suitable BEs to correct the
variation.
    <br />
    <a href="http://danioffenlab.pythonanywhere.com/"><strong>Go to website >> <strong></a>
    <br />
    <br />
    <a href="https://github.com/github_username/repo">View Demo</a>
    ·
    <a href="https://github.com/github_username/repo/issues">Report Bug</a>
    ·
    <a href="https://github.com/github_username/repo/issues">Request Feature</a>
  </p>




<!-- TABLE OF CONTENTS -->
## Table of Contents

* [About](#about-the-project)
  * [Built With](#built-with)
* [Methods to enter Dara](#getting-started)
  * [Prerequisites](#prerequisites)
  * [Installation](#installation)
* [Usage](#usage)
* [Contact](#contact)




<!-- ABOUT THE PROJECT -->
## About The Project

[![Product Name Screen Shot][product-screenshot]](https://example.com)

BE-FF is a web-based tool that receives SNV data and matches suitable BEs to correct the
variation. The code for the online tool is available here, as well as batch-mode code that may be used
to generate results for a large number of SNV's.
The main 3 scripts are:
1. <strong> BEsingle.py: </strong>
This is the main script of the website. It receives a single SNV, via 3 possible methods: <br> 
* <i> Manually entered by user:</i> <br> Here you must enter a 51-nt long DNA sequence, 25 nt upstream to the mutation and 25 nt downstream, as well as the variation and the reading frame. 


* <i> Fetched from given rsID: </i> Enters a known rsID. You will then be presented with a table containing all the possible variations. Once you selects one of the options, it will automatically be inserted into the format in (a). 


* <i> Fetched by genomic coordinates: </i> You may select a genome, chromosome number, mutation position, variation, and reading frame, which will be inserted into the format in (a).
**Note: Only one of these methods is required each time. **

You may also use the 'Advanced Options' button and create you own BE. This base editor will appear as
'User customized BE' in the final results

Press the <u> submit </u> button to show the results. You will be presented with two tables:
1. Table of  base editors that will correct this variance and no other bases around it.
2. Table of base editors that will correctly edit the SNP. However it may also edit other flanking nucleotides.
Such changes are synonymous substitutions and the resulted amino acid sequence will probably be as good as the reference sequence.





<!-- GETTING STARTED -->
## Getting Started

To get a local copy up and running follow these simple steps.


### Installation
 
1. Clone the repo
```sh
git clone https://github.com/RoyRabinowitz/BE-FF
```
2. Install NPM packages
```sh
biopython install
```



<!-- USAGE EXAMPLES -->
## Usage




<!-- CONTACT -->
## Contact

[Roy Rabinowitz](rabinowitz.roy@gmail.com) - email

[Shiri Almog](shirialmog1@gmail.com)  -email

Project Link: [https://github.com/RoyRabinowitz/BE-FF](https://github.com/RoyRabinowitz/BE-FF)



<!-- ACKNOWLEDGEMENTS -->
## Acknowledgements

* []()
* []()
* []()





<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->

[product-screenshot]: images/screenshot.png
