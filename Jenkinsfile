def modules = [
	"base",
	"debugging",
	"fieldlines",
	"galaxy",
	"globebrowsing",
	"iswa",
	"kameleon",
	"kameleonvolume",
	"multiresvolume",
	"newhorizons",
	"onscreengui",
	"space",
	"toyvolume",
	"volume"
];

def flags = "-DGHOUL_USE_DEVIL=OFF "

for (module in modules) {
	flags += "-DOPENSPACE_OPENSPACE_MODULE_" + module.toUpperCase() + "=ON "
}

echo flags

stage('Build') {
	parallel linux: {
		node('linux') {
			checkout scm
			sh 'git submodule update --init --recursive'
			sh '''
				mkdir -p build
				cd build 
				cmake .. ''' +
				flags + ''' ..
			make
			'''
		}
	},
	windows: {
		node('windows') {
			checkout scm
			bat '''
				git submodule update --init --recursive
				if not exist "build" mkdir "build"
				cd build
				cmake -G "Visual Studio 14 2015 Win64" .. ''' +
				flags + ''' ..
				msbuild.exe OpenSpace.sln /m:8 /p:Configuration=Debug
			'''
		}
	},
	osx: {
		node('osx') {
			checkout scm
			sh 'git submodule update --init --recursive'
			sh '''
				export PATH=${PATH}:/usr/local/bin:/Applications/CMake.app/Contents/bin
				export CMAKE_BUILD_TOOL=/Applications/CMake.app/Contents/bin/CMake
				srcDir=$PWD
				if [ ! -d ${srcDir} ]; then
				  mkdir ${srcDir}
				fi
				if [ ! -d ${srcDir}/build ]; then
				  mkdir ${srcDir}/build
				fi
				cd ${srcDir}/build
				/Applications/CMake.app/Contents/bin/cmake -G Xcode -D NASM=/usr/local/Cellar/nasm/2.11.08/bin/nasm ${srcDir} .. ''' +
				flags + '''
				xcodebuild
				'''
		}
	}
}	