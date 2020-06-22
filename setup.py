from setuptools import setup, find_packages
from glob import glob

setup(
    name='annotub',
    license='LICENSE',
    packages=find_packages('.'),
    package_dir={'': '.'},
    zip_safe=False,
    entry_points=dict(console_scripts=['annotub = annotub.annotub:main', ]),
    package_data={'annotub': ['*.sh']},
    data_files=[
        ('annotub/resources', glob("resources/*")),
    ],
    install_requires=['biopython', ],
    test_suite='nose.collector',
    tests_require=['nose', ],
    version='1.0'
)
