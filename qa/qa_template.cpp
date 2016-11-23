#include "qa_template.h"


void Test::set_directory(const std::string &caseFolder, const std::string &fileName)
{
	_inputDir = caseFolder + fileName + "/input/";
	_outputDir = caseFolder + fileName + "/output/";
}

std::string Test::input_file(const std::string& file)
{
	return _inputDir + file;
}
std::string Test::output_file(const std::string& file)
{
	return _outputDir + file;
}