from setuptools import setup, find_packages

setup(
    name='annotub',
    license='LICENSE',
    packages=find_packages('.'),
    package_dir={'': '.'},
    zip_safe=False,
    entry_points=dict(console_scripts=['annotub = annotub.annotub:main', ]),
    install_requires=['argparse', 'biopython', ],
    test_suite='nose.collector',
    tests_require=['nose', ],
    version='1.0'
)
