import setuptools

if __name__ == "__main__":
    setuptools.setup(
        name='mcu',
        version="1.0.1",
        description='Modeling and Crystallographic Utilities',
        project_description='A package for periodic wavefunction and crystallography analysis. mcu is designed to support large scale analysis and topological descriptions for periodic wavefunction.',
        author='Hung Pham',
        author_email='pqh3.14@gmail.com',
        url="https://github.com/hungpham2017/mcu",
        license='Apache 2.0',
        packages=setuptools.find_packages(),
        install_requires=[  
            'numpy>=1.15.2',
            'scipy>=1.1.0',
            'matplotlib>=3.0.1', 
            'spglib>=1.15.1',
        ],

        tests_require=[
            'pytest',
            'pytest-cov',
            'pytest-pep8',
            'tox',
        ],


        classifiers=[
            'Development Status :: 4 - Beta',
            'Intended Audience :: Science/Research',
            'Programming Language :: Python :: 3',
        ],
        zip_safe=True,
    )
