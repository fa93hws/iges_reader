#pragma once
#include <iostream>

class Test
{
public:
	void set_directory(const std::string &caseFolder, const std::string &fileName);

protected:
	virtual std::string input_file(const std::string& file);
	virtual std::string output_file(const std::string& file);
protected:
	std::string _inputDir;
	std::string _outputDir;
};