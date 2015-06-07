require "UTILS.rb"
def pearsonCorrelateFiles(f1,f2)
	tempFP = tempFile("R");
	tempFile = File.open(tempFP, "w");
	tempFile.print("x=read.csv(\"#{f1}\", header=FALSE);\n");
	tempFile.print("y=read.csv(\"#{f2}\", header=FALSE);\n");
	tempFile.print("result = cor.test(x=x[,1],y=y[,1], method=\"pearson\");\n");
	tempFile.print("result[['estimate']]\n");
	tempFile.print("result[['p.value']]\n");
	tempFile.close();
	output = `Rscript #{tempFP}`;
	#p(output);
	cleanFile(tempFP);
	data = output.split("\n");
	r=data[1].to_f();
	p=data[2].split("\s")[1].to_f();
	return [r,p];
end

def callRFunction(func, params);
	joinParams = "";
	if params.is_a?(Hash)
		params.keys.each(){|k|
			joinParams = joinParams+k+"=";
			if params[k].is_a?(Array)
				joinParams = joinParams+"c("+params[k].join(",")+"),";
			else
				joinParams = joinParams+params[k].to_s()+",";
			end
		}
	elsif params.is_a?(Array)
		params.each{|e|
			if e.is_a?(Hash)
				e.keys.each(){|k|
					joinParams = joinParams+k+"=";
					if e[k].is_a?(Array)
						joinParams = joinParams+"c("+e[k].join(",")+"),";
					else
						joinParams = joinParams+e[k].to_s()+",";
					end
				}
			elsif e.is_a?(Array)
				joinParams = joinParams+"c("+e.join(",")+"),";
			else
				joinParams = joinParams+e.to_s()+",";
			end
		}
	else
		joinParams = joinParams+params.to_s();
	end
	if joinParams!=""
		joinParams = joinParams[0,joinParams.length()-1];
	end
	tempFP = tempFile("R");
	`echo \"#{func}(#{joinParams})\" > #{tempFP}`
	output = `Rscript #{tempFP}`;
	cleanFile(tempFP);
	allOutput = output.split("\n");
	parsedOut = [];
	allOutput.each{|line|
		thisLine = line.split("\s");
		if thisLine[0]=~/\[\d+\]/
			parsedOut.push(thisLine[1, thisLine.length()-1]);
		end
	}
	return parsedOut;
end
