def parseSwitches(args, possible)
	switches = Hash.new();
	switch=nil;
	data=nil;
	args.each{|arg|
		if arg[0,1]=="-" && !possible[arg[1,arg.length-1]].nil?
			switch=arg[1,arg.length-1];
			data=[];
		else
			data=arg;
		end
		
		if switch.nil?
			raise("Invalid switch argument #{arg}");
		end
		if switches[switch].nil? || (switches[switch].kind_of?(Array) && switches[switch].length()==0)
			switches[switch]=data;
		else
			if !switches[switch].kind_of?(Array)
				switches[switch]=[switches[switch]];
			end
			switches[switch].push(data);
		end
	}
	#if switch!=nil
	#	switches[switch]=data;
	#end
	return switches;
end

def getDefaults(switches, defaults)
	defaults.keys().each(){|k|
		switches[k]=defaults[k] if switches[k].nil?;
	}
	return switches;
end
