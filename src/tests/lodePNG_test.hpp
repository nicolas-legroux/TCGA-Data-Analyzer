/*
 * lodePPNG_test.hpp
 *
 *  Created on: May 7, 2015
 *      Author: nicolas
 */

#ifndef SRC_TESTS_LODEPNG_TEST_HPP_
#define SRC_TESTS_LODEPNG_TEST_HPP_

#include "../lodePNG/lodepng.h"

void encodeOneStep(const char* filename, std::vector<unsigned char>& image, unsigned width, unsigned height)
{
  //Encode the image
  unsigned error = lodepng::encode(filename, image, width, height);

  //if there's an error, display it
  if(error) std::cout << "encoder error " << error << ": "<< lodepng_error_text(error) << std::endl;
}


#endif /* SRC_TESTS_LODEPNG_TEST_HPP_ */
