from setuptools import setup, find_packages
from glob import glob

setup(
    name='hybran',
    license='LICENSE',
    packages=find_packages('.'),
    package_dir={'': '.'},
    zip_safe=False,
    entry_points=dict(console_scripts=['hybran = hybran.hybran:main', ]),
    package_data={'hybran': ['*.sh']},
    data_files=[
        ('hybran/resources', glob("resources/*")),
    ],
    install_requires=['biopython', ],
    test_suite='nose.collector',
    tests_require=['nose', ],
    version='1.1.1'
)
