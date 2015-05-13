project(CancerClustering)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

file(
	GLOB_RECURSE
	src/*
)

add_executable(
	CancerClustering
	${source_files}
)
