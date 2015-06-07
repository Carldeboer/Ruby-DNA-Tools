require "fileutils"
def tempFile(ext)
	if ENV["TMPDIR"]==nil;
		tempDir=".";
	else
		tempDir = ENV["TMPDIR"];
	end

	tempFile=tempDir+"/"+"temp_"+rand(100000).to_s+"."+ext;
	while File.exists?(tempFile)
		tempFile=tempDir+"/"+"temp_"+rand(100000).to_s+"."+ext;
	end
	system("cat /dev/null >"+tempFile);
	return tempFile
end

def tempDir(prefix)
	tempDir=prefix+"temp_"+rand(100000).to_s;
	while File.directory?(tempDir)
		tempDir=prefix+"temp_"+rand(100000).to_s;
	end
	return tempDir
end

def cleanFile(filePath)
	if File.exists?(filePath)
		File.delete(filePath);
	end
end
def cleanDir(dirPath)
	if File.directory?(dirPath) && dirPath!="." && dirPath!=".."
		FileUtils.rm_rf(dirPath);
	end
end


def fileDiff(fp1, fp2)
	if !File.exists?(fp1) || !File.exists?(fp2)
		return true;
	end
	return !FileUtils.compare_file(fp1, fp2);
end


def html2xml(html)
	tempFile = tempFile("html");
	htmlFile = File.open(tempFile,"w");
	htmlFile.print(html);
	htmlFile.close();
	#xml = `tidy -asxml #{tempFile} `;
	xml = `tidy -asxml #{tempFile} 2>/dev/null`;
	return xml;
	cleanFile(tempFile);
end
