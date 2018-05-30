# DomainSeq

# Overview

## Requirements
This project has only been tested on Unix/Linux (OSX,Ubuntu,etc). The project is not platform specific, but installation instructions only apply to these systems. Additionally, all code is formated for Python 3.X.

## Installation
We recommend using this package through a [virtual environment](https://docs.python.org/3/library/venv.html). This creates an "instance" of Python that is isolated from the rest of your system. To do this, open terminal:

1. Clone this package
1. Install the virtual environment package to your native Python:

   > python3 -m pip install --user virtualenv

2. Use this package to make a virtual environment of Python 3:

   > python3 -m virtualenv domainseq-venv
   
   Note: This creates a folder in your current directory called "domainseq-venv" which contains your virtual environment. Call it whatever you want!
   
3. Activate your virtual environment by declaring it as your source for Python:

   > source domainseq-venv/bin/activate

   Note: This command causes the virtual environment (your "domainseq-venv" folder) to be prioritized as the version of Python your session will use.
   
4. If you are not already there, navigate to the home of this folder
