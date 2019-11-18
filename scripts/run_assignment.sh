set -e

if [ -z "$1" ]
then
	echo >&2 "No assignment given.";
	exit 1;
fi

assignment=$1
test=$2

cmake --build ./build --config Release --target opengl_viewer
cmake --build ./build --config Release --target $assignment

./build/proj/$assignment/$assignment -test $test

