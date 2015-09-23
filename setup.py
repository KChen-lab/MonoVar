import os
import sys
import subprocess
try:
    from setuptools import setup, Extension
    from setuptools.command.install import install
    from setuptools.command.develop import develop
    from distutils.command.build import build
    havesetuptools = True
except ImportError:
    from distutils.core import setup, Extension
    from distutils.command.install import install
    from distutils.command.build import build
    havesetuptools = False

BASEPATH=os.path.dirname(os.path.abspath(__file__))

class MonoVarBuild(build):

	def run(self):
		# run original build code
		build.run(self)

		# build samtools
		build_path = os.path.abspath(self.build_temp)
		cmd = ['make', '-C', 'external/samtools']

		def compile():
			subprocess.check_call(cmd)

		self.execute(compile, [], 'Compile samtools')

class MonoVarInstall(install):

	def run(self):
		install.run(self)
		import shutil
		shutil.copy2('external/samtools/samtools',
                     os.path.join(self.install_lib, 'monovar_src'))

cmdclass = {'build': MonoVarBuild, 'install': MonoVarInstall,}

def main():
	if float(sys.version[:3])<2.6 or float(sys.version[:3])>=2.8:
        	sys.stderr.write("CRITICAL: Python version must be 2.6 or 2.7!\n")
        	sys.exit(1)

	setup(
	name = "MonoVar",
	author = "Hamim Zafar",
	cmdclass = cmdclass,
	)

if __name__ == '__main__':
    main()
