def smartGZIn(theFP,mode)
	if theFP[-3,3].upcase()==".GZ"
		return inGZ(theFP);
	else
		File.open(theFP,"r");
	end
end	
def smartGZOut(theFP,mode)
	if theFP[-3,3].upcase()==".GZ"
		return outGZ(theFP,mode);
	else
		File.open(theFP,mode);
	end
end	

def smartGZForeach(theFP)
	if block_given?
		if theFP[-3,3].upcase()==".GZ"
			inGZ(theFP).each{|line|
				yield(line);
			}
		elsif theFP[-4,4].upcase()==".ZIP"
			inZIP(theFP).each{|line|
				yield(line);
			}
			#return readGZ(theFP).split("\n").map();
		else
			File.foreach(theFP){|line|
				yield(line);
			}
		end
	else
		if theFP[-3,3].upcase()==".GZ"
			return inGZ(theFP);
		elsif theFP[-4,4].upcase()==".ZIP"
			return inZIP(theFP);
		else
			return File.foreach(theFP);
		end
	end
end

def inZIP(inFP)
	return IO.popen("unzip -p #{inFP}");
end

def inGZ(inFP)
	return IO.popen("gunzip -c #{inFP}");
end
def outGZ(outFP,mode="w")
	return IO.popen("gzip -c >#{outFP}",mode);
end

def readGZ(inFP)
	return `gunzip -c #{inFP}`
end
