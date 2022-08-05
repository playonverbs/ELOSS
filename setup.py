from distutils.core import setup
setup(
    name = 'ELOSS',
    packages = ['ELOSS'],
    version = '0.1.0',
    license = 'None',
    description = 'A library for modelling energy loss in argon',
    author = 'David Caratelli',
    author_email = 'user@example.com',
    url = 'https://github.com/davidc1/ELOSS',
    download_url = 'https://github.com/davidc1/ELOSS/archive/refs/head/master.tar.gz',
    keywords = ['nuclear modelling', 'argon'],
    install_requires = ['numpy >= 1.15.0'],
    extras_require = {
        'tests': [
            "matplotlib",
            "pandas",
            "scipy",
        ]
    },
    classifiers = [
        'Intended Audience :: Developers',
        'Programming Language :: Python :: 2',
    ],
)
