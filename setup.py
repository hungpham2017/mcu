import setuptools

if __name__ == "__main__":
    setuptools.setup(
        name='mcu-hungpham2017',
        version="0.0.2",
        description='Modeling and Crystallographic Utilities',
        author='Hung Pham',
        author_email='pqh3.14@gmail.com',
        url="https://github.com/hungpham2017/mcu.git",
        license='Apache 2.0',
        packages=setuptools.find_packages(),
        install_requires=[
            'python>=3.*',        
            'numpy>=1.15.2',
            'scipy>=1.1.0',
            'matplotlib>=3.0.1', 
            'spglib>=1.12.0',
        ],
        extras_require={
            'docs': [
                'sphinx==1.2.3',
                'sphinxcontrib-napoleon',
                'sphinx_rtd_theme',
                'numpydoc',
            ],
            'tests': [
                'pytest',
                'pytest-cov',
                'pytest-pep8',
                'tox',
            ],
        },

        tests_require=[
            'pytest',
            'pytest-cov',
            'pytest-pep8',
            'tox',
        ],


        classifiers=[
            'Development Status :: 4 - Beta',
            'Intended Audience :: Science/Research',
            'Programming Language :: Python :: 3.7.6',
        ],
        zip_safe=True,
    )
