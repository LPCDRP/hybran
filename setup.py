from setuptools import setup, find_packages

setup(
    name='hybran',
    license='LICENSE',
    packages=find_packages('src'),
    package_dir={'': 'src'},
    zip_safe=False,
    entry_points=dict(console_scripts=['hybran = hybran.hybran:main', ]),
    package_data={'hybran': [
        '*.sh',
    ]},
    install_requires=['biopython', ],
    test_suite='nose.collector',
    tests_require=['nose', ],
    version='1.5.2'
)
