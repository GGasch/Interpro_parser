from setuptools import setup

setup(
    name='interpro_parser',
    version='1.0.0',    
    description='A code to extract domain sequences based on interpro data (table file and fasta file)',
    url='https://github.com/GGasch/interpro_parser.git',
    author='Geoffroy Gaschignard',
    author_email='geoffroy.gaschignard@gmail.com',
    license='MIT Licence',
    packages=['interpro_parser'],
    install_requires=['Bio',
                      'tqdm',
                      'argparse', 
                      'pandas',                   
                      ],

    classifiers=[
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT Licence',  
        'Operating System :: POSIX :: Linux',        
        'Programming Language :: Python :: 3',
   
    ],
)

