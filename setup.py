from setuptools import setup, find_packages

setup(name='sigmd',
      version='0.1',
      description='CMDS signal analysis',
      url='',
      author='S. Palato',
      author_email='',
      license='',
      packages=find_packages("src"),
      package_dir={"": "src"},
      zip_safe=False)