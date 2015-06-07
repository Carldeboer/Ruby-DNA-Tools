
#ref\tsyn\tsyn\tsyn...
def loadDictionary(filePath)
	thisHash = Hash.new()
	File.foreach(filePath){|line|
		line.chomp!
		splitLine = line.split("\t")
		reference = splitLine.shift().upcase()
		thisHash[reference]=reference
		splitLine.each{|curKey|
			curKey.upcase!()
			if thisHash[curKey]!=nil
				if thisHash[curKey]!=reference
					print("Conflicting entry for "+curKey+"->"+thisHash[curKey]+" and "+curKey+"->"+reference+"\n")
				end
			end
			thisHash[curKey]=reference
		}
	}
	return thisHash;
end

#key->[data]
def loadHashFile(filePath)
	thisHash = Hash.new()
	File.foreach(filePath){|line|
		line.chomp!
		splitLine = line.split("\t")
		key = splitLine.shift()
		thisHash[key]=splitLine
	}
	return thisHash;
end
