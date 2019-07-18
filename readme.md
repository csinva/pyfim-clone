# pyfim clone

Clone of pyfim making it installable as a dependency. Copied from http://www.borgelt.net/pyfim.html - I did not write any of this code.

# installation
`pip install git+https://github.com/csinva/pyfim-clone`

Now we can import pyfim (from python):

`from fim import fpgrowth`

# include in your project

If you would like to include pyfim as a dependency for your project, simply add this to your `setup.py` file:
```python    
setup(
    ...
    install_requires=[
        'fim',
    ],
    dependency_links=[
        'https://github.com/csinva/pyfim-clone/tarball/master#egg=fim-6.28'
    ],
    ...
)
```
