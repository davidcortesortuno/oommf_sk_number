import setuptools

with open('README.md', encoding='utf-8') as f:
    long_description = f.read()

setuptools.setup(
    name='oommf_sk_number',
    version='0.1',
    description=('Scripts for the calculation of sk number from '
                 'OOMMF micromagnetic code files'),
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='D. Cortes',
    author_email='d.cortes@soton.ac.uk',
    packages=setuptools.find_packages(),
    install_requires=['matplotlib',
                      'numpy'],
    classifiers=['License :: BSD2 License',
                 'Programming Language :: Python :: 3 :: Only',
                 ]
)
