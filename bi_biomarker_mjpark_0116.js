"use strict"

/***********************************************
	patch log
 ***********************************************
 2017.03.02	kyuwan	remove //_log_debug
 2017.03.03 kyuwan	add	function				_make_cladogram_data()
 2017.03.08 kyuwan	modified for groups(n>2)	_make_cladogram_data()
 2017.05.23 jeongm	add function to select 500 features
 2017.06.20 minje	modify function to select 500 features logic
 2017.07.11 minje	modify lefse main logic
 						- remove target_num
 						- remove bootstrapping option (n_boots=1,fract_sample=1)
 						- remove line 929 	normal()
 						- change implementation gamma/log-gamma functions, using mathjs
 						- change getNormalZ --> p_to_z using mathjs
2017.07.13 minje	modify
 						- revive target_num (for speed issue)
 						- fix sort error changed in 0620	
2017.07.14 minje	modify
 						- revive line 929 (for singular matrix issue)
 						 	- change stdev 0.05 --> 0.000000005
 						- #feature >= 500 , prune to 499
2017.11.06 minje 	modify
						- when data goes to singular throw error 

2017.11.29 minje	rework
						- wilcoxon,kruskalwallis,xor
					add functon
						- wilcoxon_core,kruskal_core,xor_core
					modify lda module
						- math --> numeric

2017.11.30 minje	modify 
						- error handling when NaN/Infinity/undefined/empty comes....

2017.12.04 minje	fix
						- kruskal_core logic if input elements are all same...
					add
						- log for some steps...(for debug...)

2017.12.12 minje	fix
						- console.log / console.error --> _log_debug
						- insert input parameter "server_log_debug"

2017.12.20 minje	fix
						- xor_core logic (minor)
					add function
						- permanova //actually this function should be on "beta-diversity...."
							-permanova
							-permanova_core
							-getF
							-get_ssw

2017.12.26 minje	modify
						-permanova return value 
							insert:  basic statistics (stat_for_boxplot)
							remove:  distance matrix itself...
					add function
						-stat_for_boxplot --> returns q1 q2 q3 avg stdev variance min max

2018.01.09 minje	fix
						-permanova output
2018.01.18 minje	fix
						-permanova boxplot data (q1,q3) --> mathjs interpolates quantile data....

--> change inverse method as pseudo-inverse later on....
***********************************************/

var math = require('mathjs');
var numeric = require('numeric');

var _log_debug = function(msg) {console.log('>> ' + msg) };

function wilcoxon(input,server_log_debug){
	if (server_log_debug != undefined && server_log_debug != null) {
		_log_debug = server_log_debug;
	}

	if(input.group.unique().length !== 2){
		_log_debug("#group must be 2 in wilcoxon");
	}else{
		let in_sample = json_to_array(input);
		_log_debug("[info] json --> array complete :"+math.size(in_sample.profile))
		const root_count = in_sample.profile.map(function(v,i){return v[0]});
		let res = []
		_log_debug("[info] calculate wilcoxon for "+in_sample.profile[0].length+" nodes")
		for(let i=1;i<in_sample.profile[0].length;i++){ // as i=0 --> root, we skip root
			let in_array = in_sample.profile.map(function(val,idx){return val[i]/root_count[idx]});
			let in_array_split = in_array.chunk(wordcount(in_sample.group).count);
			res.push({key : in_sample.keys[i], p: wilcoxon_core.apply(this,in_array_split).pvalue});
		}
		getQ(res);
		_log_debug("[info] calculation complete")

		let res_005 = res.filter(function(v){return parseFloat(v.p)<0.05});
		let Result=[];
		for(let i=0;i<res_005.length;i++){
			let block = res_005[i];
			if(block.key <= 0) continue;
			let ratio =  in_sample.profile.map(function(val,idx){return val[in_sample.keys.indexOf(block.key)]/root_count[idx]}); //long...
			Result.push({tn : in_sample.keys[in_sample.keys.indexOf(block.key)], hierarchy : getParents(block.key,in_sample), nm : in_sample.name[in_sample.keys.indexOf(block.key)],
						lvl : Number(in_sample.level[in_sample.keys.indexOf(block.key)]), pValue: [block.p,block.q],
						ratio : ratio,average : ratio.chunk(wordcount(in_sample.group).count).map(function(v){return v.reduce(function(a,b){return a+b})/v.length }) });
		}
		let Result_wilcoxon = Result.filter(function(v){return v.ratio.reduce(function(a,b){return a+b}) > 0.0001});
		_log_debug("[info] result --> json")
		return {wilcoxonResult : Result_wilcoxon, group_names : input.groupName};
	}
}

function kruskal(input,server_log_debug){
	if (server_log_debug != undefined && server_log_debug != null) {
		_log_debug = server_log_debug;
	}

	if(input.group.unique().length <2){
		_log_debug('#group must be equal or greater than 2 in kruskal');
	}else{
		let in_sample = json_to_array(input);
		_log_debug("[info] json --> array complete :"+math.size(in_sample.profile))
		let root_count = in_sample.profile.map(function(v,i){return v[0]});
		let res = []
		_log_debug("[info] calculate kruskal for "+in_sample.profile[0].length+" nodes")
		for(let i=1;i<in_sample.profile[0].length;i++){	//as i=0 --> root, we skip root
			let in_array = in_sample.profile.map(function(val,idx){return val[i]/root_count[idx]});
			let in_array_split = in_array.chunk(wordcount(in_sample.group).count);
			res.push({key : in_sample.keys[i], p: kruskal_core.apply(this,in_array_split).pvalue});
		}
		getQ(res);
		_log_debug("[info] calculation complete")
		
		let res_005 = res.filter(function(v){return v.p<0.05});
		let Result=[];
		for(let i=0;i<res_005.length;i++){
			let block = res_005[i];
			let ratio =  in_sample.profile.map(function(val,idx){return 100*val[in_sample.keys.indexOf(block.key)]/root_count[idx]}); //long...
			Result.push({tn : block.key, hierarchy : getParents(block.key,in_sample), nm : in_sample.name[in_sample.keys.indexOf(block.key)],
						lvl : in_sample.level[in_sample.keys.indexOf(block.key)], pValue: [block.p,block.q],
						ratio : ratio,average : ratio.chunk(wordcount(in_sample.group).count).map(function(v){return v.reduce(function(a,b){return a+b})/v.length}) })
		}
		let Result_kruskal = Result.filter(function(v){return v.ratio.reduce(function(a,b){return a+b}) > 0.0001})
		_log_debug("[info] result --> json")
		return {kruskalResult : Result_kruskal, group_names : input.groupName};
	}
}

function xor(input,server_log_debug){
	if (server_log_debug != undefined && server_log_debug != null) {
		_log_debug = server_log_debug;
	}

	if(input.group.unique().length <2){
		_log_debug('#group must be equal or greater than 2 in xor');
	}else{
		var in_sample = json_to_array(input);
		_log_debug("[info] json --> array complete :"+math.size(in_sample.profile))
		const root_count = in_sample.profile.map(function(v,i){return v[0]});
		let res = []
		_log_debug("[info] calculate xor for "+in_sample.profile[0].length+" nodes")
		for(let i=0;i<in_sample.profile[0].length;i++){
			let in_array = in_sample.profile.map(function(val,idx){return val[i]});
			let in_array_split = in_array.chunk(wordcount(in_sample.group).count);
			res.push({key : in_sample.keys[i], xor: xor_core.apply(this,in_array_split), arr: in_array});
		}
		_log_debug("[info] calculation complete")
		let res_true = res.filter(function(v){return v.xor !== false});
		let Result = [];
		for(let i=0;i<res_true.length;i++){
			let block = res_true[i]
			let ratio =  in_sample.profile.map(function(val,idx){return 100*val[in_sample.keys.indexOf(block.key)]/root_count[idx]}); //long...
			Result.push({tn : block.key, hierarchy : getParents(block.key,in_sample), nm : in_sample.name[in_sample.keys.indexOf(block.key)],
						lvl : in_sample.level[in_sample.keys.indexOf(block.key)], ratio : ratio, existence: block.xor.existence, 
						specific : block.xor.specific
					})
		}
		_log_debug("[info] result --> json")
		return {xorResult: Result,group_names : input.groupName};
	}
}

function wilcoxon_core(input){ //return p value...
	const args =  Array.from(arguments);
	if(args.length !== 2){
		_log_debug("[E] Wilcoxon need 2 input");
	}else{
		//in scipy, it only use this, regardless of tie...
		const arr = args.reduce(function(p,c){return p.concat(c)});
		if(Object.keys(arr).length<arr.length){ //there is empty element
			_log_debug("[E] your input contains empty element")
		}else{
			if(arr.every(function(elem){return isFinite(elem)})==false){
				_log_debug("[E] your input contains NaN/Infinity/undefined")
			}else if(arr.some(function(elem){return elem == null})){
				_log_debug("[E] your input contains null value")
			}else{
					let arr_sort = arr.map(Math.abs).slice().sort(function(a,b){return a-b;});
					let arr_rank = arr.map(Math.abs).slice().map(function(v){return arr_sort.indexOf(v)+1;});
					let rank_mod = rankTie(arr_rank.slice());
					let rank_split = rank_mod.slice().chunk(args.map(function(v){return v.length}));
					let s = rank_split[0].reduce(function(a,b){return a+b});
					let expected = rank_split[0].length * (rank_mod.length+1)/2
					let z = (s-expected) / math.sqrt(rank_split[0].length*rank_split[1].length*(rank_mod.length+1)/12)
					let p = z_to_p(math.abs(z),"two-tailed")
					return {statistics : z, pvalue: p}
			}
		}
	
	}
}

function kruskal_core(input){
	const args =  Array.from(arguments);
	if(args.length<2){
		_log_debug("[E] Kruskal need more or equal than 2 input");
	}else{
		const arr = args.reduce(function(p,c){return p.concat(c)});
		if(Object.keys(arr).length<arr.length){ //there is empty element
			_log_debug("[E] your input contains empty element")
		}else{
			if(arr.every(function(elem){return isFinite(elem)})==false){
				_log_debug("[E] your input contains NaN/Infinity/undefined")
			}else if(arr.some(function(elem){return elem == null})){
				_log_debug("[E] your input contains null value")
			}else{
				let arr_sort = arr.map(Math.abs).slice().sort(function(a,b){return a-b;});
				let arr_rank = arr.map(Math.abs).slice().map(function(v){return arr_sort.indexOf(v)+1;});
				let rank_mod = rankTie(arr_rank.slice());
				let rank_split = rank_mod.slice().chunk(args.map(function(v){return v.length}));
				let r_bar = 0.5*(arr.length+1);
				let g_rank_sum = rank_split.map(function(v){return v.reduce(function(p,c){return p+c})});

				if(wordcount(rank_mod).count.every(function(elem){return elem == 1})){ //no duplicate
					let H_noTie = (12/(arr.length*(arr.length+1))) * g_rank_sum.map(function(v,i){return v*v/args[i].length}).reduce(function(p,c){return p+c}) - 3*(arr.length+1);
					let p = getChiCdf(H_noTie,args.length-1);
					if(isNaN(p)) p=1;
					return {statistics: H_noTie,pvalue:p};

				}else{
					let H_noTie = (12/(arr.length*(arr.length+1))) * g_rank_sum.map(function(v,i){return v*v/args[i].length}).reduce(function(p,c){return p+c}) - 3*(arr.length+1);
					let denominator = 1 - (wordcount(rank_mod).count.filter(function(elem){return elem >1}).map(function(v){return (Math.pow(v,3)-v)}).reduce(function(p,c){return p+c})/(Math.pow(arr.length,3)-arr.length));
					
					if(denominator==0){
						return {statistics: 0,pvalue: 1};
					}else{
						let H_tie = H_noTie / denominator;
						let p = getChiCdf(H_tie,args.slice().length-1);
						if(isNaN(p)) p=1;
						return {statistics: H_tie,pvalue:p};
					}
				}
			}
		}
	}
}

function xor_core(input){
	const args =  Array.from(arguments);
	if(args.length<2){
		_log_debug("[E] xor need more or equal than 2 input");
	}else{
		const arr = args.reduce(function(p,c){return p.concat(c)});
		if(Object.keys(arr).length<arr.length){ //there is empty element
			_log_debug("[E] your input contains empty element")
		}else{
			if(arr.every(function(elem){return isFinite(elem)})==false){
				_log_debug("[E] your input contains NaN/Infinity/undefined")
			}else if(arr.some(function(elem){return elem == null})){
				_log_debug("[E] your input contains null value")
			}else{
				//check group --> 0/1
				let args_zero = args.map(function(v){return v.every(function(elem){return elem==0})})    //all zero
				let args_nzero = args.map(function(v){return v.every(function(elem){return elem!==0})})   //all non zero
				let arg_flag = numeric.add(args_nzero,args_zero);
				
				if(arg_flag.every(function(elem){return elem == 1})){
					let flag = numeric.add(args_nzero,numeric.mul(0,args_zero)); //non zero --> 1, zero --> 0
					if( args.length > 2 && wordcount(flag).count.filter(function(elem){return elem==1}).length==1 ){
						let spec= flag.indexOf(wordcount(flag).word[wordcount(flag).count.indexOf(1)]) 
						let res = arr.slice().map(function(v){if(v == 0){return false}else{return true}})
						return {existence : res, specific : spec}

						}else if(args.length == 2  && wordcount(flag).word.length ==2){
							let spec= flag.indexOf(1) //when comparing two groups, specific group becomes "1"
							const arr = args.reduce(function(p,c){return p.concat(c)});
							let res = arr.slice().map(function(v){if(v == 0){return false}else{return true}})
							return {existence : res, specific : spec}
						}else{
							return false
						}
				}else{
					return false
				}
			}
		// it seems that [1 0 1 0] returns true, but former xor returns false to specify return variable "specific"
		}
	}
}



function json_to_array(input){
	var group = input.group;
	var json = input.samples;
	var key = new Array();
	for(let i=0;i<json.length;i++){
		key.push(json[i].map(function(v){return v["tn"]}));
	}
	var keys = ([].concat.apply([],key)).unique();
	var cntTable = math.zeros([json.length,keys.length]);
	var parent = [];
	var name = [];
	var level = [];

	for(let i=0;i<json.length;i++){
		for(let j=0;j<json[i].length;j++){
			let tmp = json[i][j]["cnt"][0];
			cntTable[i][Number(keys.indexOf(json[i][j]["tn"]))] = tmp;
			parent[Number(keys.indexOf(json[i][j]["tn"]))] = json[i][j]["p"];
			name[Number(keys.indexOf(json[i][j]["tn"]))] = json[i][j]["nm"];
			level[Number(keys.indexOf(json[i][j]["tn"]))] = json[i][j]["lvl"];
		}
	}
	var res = {group : group, keys: keys, profile : cntTable, parent : parent, name : name, level : level};
	return res;
}

function getParents(index,table){
	let idx = index;
	let nestName=[];
	if(table.parent[table.keys.indexOf(idx)] !== -1){
		while(table.parent[table.keys.indexOf(idx)] !== -1){
			nestName.push(table.name[table.keys.indexOf(idx)]);
			idx = table.parent[table.keys.indexOf(idx)];
		}
		return nestName.reduceRight(function(a,b){return a+": "+b;});
	}else{
		return nestName.push('root');
	}
}

function groupFilter(input,groupData,Gvalue){
	var Gfilter = input.map(function(value,index){if(groupData[index]==Gvalue) return value});
	return Gfilter.filter(function(element){return element !== undefined});
}

// _.chunk(A,xx) --> split by only even size...

function wordcount(input){
	let count =[];
	let word = input.unique();
	for(let i=0;i<word.length;i++){
		let cc = input.filter(item => item == word[i]).length; //IE don't support ES6 element...
		//let cc = input.filter(function(item,i){return item==word[i]}).length;
		count.push(cc);
	}
	return {word : word,count : count};
}

function rankTie(input){
	let count = wordcount(input);
	return input.map(function(v,i){
		return v+(count.count[count.word.indexOf(v)]-1)*0.5;
	});
}

function getQ(input){
	//using Benjamini-Hochberg method
	//let alpha = 0.05 //doesn't really matter....
	let sorted = input.slice().sort(function(a,b){return b.p-a.p});
	sorted.forEach(function(v,i){sorted[i].q = Math.abs(v.p)*sorted.length/(sorted.length-i)});
	sorted.map(function(v){return (v>=1) ? 1 : v});
	return sorted.slice().sort(function(a,b){return a.key-b.key});
}

function normal(mu, sigma, nsamples){

	if(!nsamples) nsamples = 1
	if(!sigma) sigma = 1
	if(!mu) mu=0

	var run_total = 0
	for(var i=0 ; i<nsamples ; i++){
		run_total += Math.random();
	}

	return sigma*(run_total - nsamples/2)/(nsamples/2) + mu
}


function z_to_p(z,option){
	// 0.5 + 0.5*erf(x-mu/(sigma*root2)) mu=0,sigma=1
	if(option=='one-tailed'){
		return 1- (0.5 + 0.5*math.erf(z/math.sqrt(2)));
	}else if(option=='two-tailed'){
		return 2*(1- (0.5 + 0.5*math.erf(z/math.sqrt(2))));
	}
}


//function related gamma function is from 
//numerical recipies in C
//http://www.it.uom.gr/teaching/linearalgebra/NumericalRecipiesInC/c6-2.pdf

function Gammacdf(x,a){
	var GI;
	if (x<=0) {
		GI=0
	} else if (x<a+1) {
		GI=gser(x,a)
	} else {
		GI=gcf(x,a)
	}
	return GI
}

function getChiCdf(x,a){
    var Z=eval(x)
    var DF=eval(a)
    var Chisqcdf=0;
	if (DF<=0) {
		alert("Degrees of freedom must be positive")
	} else {
		Chisqcdf=Gammacdf(Z/2,DF/2)
	}
	//Chisqcdf=Math.round(Chisqcdf*100000)/100000;
    return Chisqcdf;
}
/*
function getChiCdf(x,a){
    var Z=eval(x)
    var DF=eval(a)
    var Chisqcdf=0;
	if (DF<=0) {
		alert("Degrees of freedom must be positive")
	} else {
		Chisqcdf=Gammacdf(Z/2,DF/2)
	}
	Chisqcdf=Math.round(Chisqcdf*Math.pow(10,32))/Math.pow(10,32);
    return Chisqcdf;
}
*/
function beta(z,w){
	return math.gamma(z)*math.gamma(w)/math.gamma(z+w);
}

function gcf(x,a){
	const max_iter = math.pow(10,12);
	const eps = 3.0*math.pow(10,-20);
	const fp_min = 1.0*math.pow(10,-60);

	let b = x + 1.0-a;
	let c = 1.0 / fp_min;
	let d = 1.0 / b;
	let h = d;
	for(let i=1;i<=max_iter;i++){
		let an = -i*(i-a);
		b += 2.0;
		d = an*d+b;
		if(math.abs(d)<fp_min) d = fp_min;
		c=b+an/c;
		if(math.abs(c)<fp_min) c = fp_min;
		d=1.0/d;
		let del=d*c;
		h *= del;
		if(math.abs(del-1.0)<eps) break;
	}
	return math.exp(-x+a*math.log(x)-math.log(math.gamma(a)))*h;
}

function gser(x,a){
	const max_iter = math.pow(10,12);
	const eps = 3.0*math.pow(10,-20);
	const fp_min = 1.0*math.pow(10,-60);

	if(x<0.0){
		_log_debug("x must be greater than 0 in gser")
	}else if(x>0.0){
		let ap=a;
		let del = 1.0/a;
		let sum = 1.0/a;
		for(let i=1;i<=max_iter;i++){
			++ap;
			del *= x/ap;
			sum += del;
			if(math.abs(del)<math.abs(sum)*eps){
				return sum*math.exp(-x+a*math.log(x)-math.log(math.gamma(a)));
			}
		}
	}else{
		return 0;
	}
}

function diag(input){
	//diag of nxn matrix
	if(input.length==input[0].length){
		return input.map(function(v,i){return v[i]});
	}else{
		_log_debug("input matrix is mxn, it must be nxn")
	}
}

Array.prototype.unique = function(){ // --> to _.uniq in lodash
   var u = {}, a = [];
   for(var i = 0, l = this.length; i < l; ++i){
      if(u.hasOwnProperty(this[i])) {
         continue;
      }
      a.push(this[i]);
      u[this[i]] = 1;
   }
   return a;
}

Array.prototype.chunk = function(input){
	// [x,x,x,x,x] --> [x,x,x],[x,x] ! not  [x,x,x,x,x] --> [x,x,x],[x,x]
	//										[y,y,y,y,y]     [y,y,y],[y,y,] 
	let res = [];
	let start = 0;
	for(let i=0;i<input.length;i++){
		let end = start + input[i];
		res.push(this.slice(start,end));
		start = start + input[i];	
	}
	return res;		
}

// added by kyuwan (2017.03.03)
function _make_cladogram_data(samples,output) {
	var cdata = {name: 'root'};

	for (var i=0;i<samples.length;i++) {
		for (var j=0;j<samples[i].length;j++) {
			if (samples[i][j].lvl == '12' || samples[i][j].nm.endsWith('_uc')) {	
				var h = [samples[i][j].nm];		
				var np = samples[i][j].np;
				while (np > 0) {
					h.push(samples[i][np].nm);
					np = samples[i][np].np;	
				}
				
				var node = cdata;
				for (var k=(h.length - 1);k>=0;k--) {
					if (node.children == undefined) {
						var child = {name:h[k]};
						node.children = [child];
						node = child;	
						continue;
					} else {
						var bFind = false;
						for (var m=0;m<node.children.length;m++) {
							if (node.children[m].name == h[k]) {
								node = node.children[m];
								bFind = true;
								break;
							}
						}
						if (bFind) {
							continue;
						}

						var child = {name:h[k]};
						node.children.push(child);
						node = child;	
					}
				}
			}
		}
	}



	for (var i=0;i<output.lefse_result.length;i++) {
		var h = output.lefse_result[i].hierarchy.split(': ');
	
		var node = cdata;
		for (var j=0;j<h.length;j++) {
			if (node.children == undefined) {
				var child = {name:h[j]};
				node.children = [child];
				node = child;	
				continue;
			} else {
				var bFind = false;
				for (var k=0;k<node.children.length;k++) {
					if (node.children[k].name == h[j]) {
						node = node.children[k];
						bFind = true;
						break;
					}
				}
				if (bFind) {
					continue;
				}

				var child = {name:h[j]};
				node.children.push(child);
				node = child;	
			}
		}

		var maxV = math.max(output.lefse_result[i].average);
		var maxI = output.lefse_result[i].average.indexOf(maxV);
		node.radius = math.log(output.lefse_result[i].average[maxI] * 10000,10);
		node.pValue = output.lefse_result[i].pValue[0];
		node.group = maxI;
	}

	return cdata;
}



//Normal distribution and Chi-square distribution code may be changed someday...

function lefse(input,server_log_debug){
	if (server_log_debug != undefined && server_log_debug != null) {
		_log_debug = server_log_debug;
	}

	var groupInput = input.group;
	var sampleN = groupInput.length;
	var groupArr = groupInput.unique();
	var groupN = groupArr.length;

	var nByGroup = new Array(groupN);
	
	for (var groupIdx = 0; groupIdx < groupN ; groupIdx++)
	{
		nByGroup[groupIdx] = 0;
	}

	for (var i = 0; i < sampleN ; i++)
	{
		nByGroup[groupInput[i]]++;
	}
	
	var groupDoubleArr = [];
	for (var i = 0; i < groupN; i++)
	{
		groupDoubleArr[i] = new Array();
	}

	for (var i = 0; i < nByGroup.length ; i++)
	{
		var n = nByGroup[i];
		for (var j = 0; j < n ; j++)
		{
			groupDoubleArr[i].push(i);
		}
	}

// 1. Kruskal Wallis
	//console.time("kw");
	var res = kruskal(input);
	//console.timeEnd("kw");
	if(res.kruskalResult.length == 0){
		var lefse_result = {lefse_result : [], group_names : input.groupName};
		return lefse_result;
	}
	var biomarkers = res.kruskalResult;
	var groupNames = res.group_names;

	var matX = new Array(groupN);

	for (var groupIdx = 0; groupIdx < groupN; groupIdx++)
	{
		matX[groupIdx] = new Array(biomarkers.length);
	}
	
	for (var i = 0; i < biomarkers.length; i++){
		var arr = biomarkers[i].ratio.map(function(value){return value*10000;});
		var tmp_arr = [];
		for (var groupIdx = 0; groupIdx < groupN ; groupIdx++)
		{
			tmp_arr[groupIdx] = [];
		}
		for (var j = 0; j < groupInput.length ; j++){
			tmp_arr[groupInput[j]].push(arr[j]);
		}
		for (var groupIdx = 0; groupIdx < groupN ; groupIdx++)
		{
			matX[groupIdx][i] = tmp_arr[groupIdx];
		}
	}

	var mean_arr = [];
	var matX_t = [];

	for (var groupIdx = 0; groupIdx < groupN; groupIdx++)
	{
		mean_arr.push(math.mean(matX[groupIdx],1));
		matX_t.push(math.transpose(matX[groupIdx]));
	}

// 2. Wilcoxon on subclasses
	//console.time("willcoxon")
	//check if i was right,,,
	
	for(let i=0;i<biomarkers.length;i++){
		let pairwise_p = []
		for(let gIdx = 0; gIdx<groupN;gIdx++){
			let p_temp = []
			for(let gRef = gIdx+1;gRef<groupN;gRef++){
				let cIdx = matX[gIdx][i]
				let cRef = matX[gRef][i]
				p_temp.push(kruskal_core(cIdx,cRef).pvalue)
			}
			pairwise_p.push(p_temp);
		}
		//_log_debug(pairwise_p.filter(function(elem){return elem<0.05}))
		//_log_debug(pairwise_p)
	}
	
	var alpha_mtc = 0.05; // significant level
	var sig_wilc = [];

	for (var i = 0; i < biomarkers.length; i++)
	{
		
		var tot_ok = 0;
		var all_diff = [];

		for ( var groupIdx = 0; groupIdx < groupN; groupIdx++ )
		{
			for( var groupIdxRef = (groupIdx + 1); groupIdxRef < groupN; groupIdxRef++ )
			{
				var first = true;
				var ok = 0;
				var br = false;

				var cl1 = matX[groupIdx][i];
				var cl2 = matX[groupIdxRef][i];
				var med_comp = false;
				if(cl1.length < 10 | cl2.length.length < 10) med_comp = true;
				var sx = math.median(cl1);
				var sy = math.median(cl2);
				if (cl1[0] == cl2[0] & set(cl1).length == 1 & set(cl2).length == 1)
				{
					var tres = false;
					first = true;
				} else if (!med_comp)
				{
					var tmp_groupInput = groupDoubleArr[groupIdx].concat(groupDoubleArr[groupIdxRef]);
					//var pv = kruskalWallis_lefse(cl1,cl2,tmp_groupInput);
					var pv = kruskal_core(cl1,cl2).pvalue;
					var tres = pv < alpha_mtc * 2.0;
				}
				if (first)
				{
					first = false;
					if (med_comp | tres)
					{
						var dir_cmp = sx < sy;
					} else
					{
						br = true;
					}
				} else if (med_comp)
				{
					if ((sx < sy) != dir_cmp | sx == sy) br = true;
				} else if (!tres | (sx < sy) != dir_cmp | sx == sy) br = true;
				if (br) continue;
				ok++;

				var diff = (ok == 1);
				if(diff){
					tot_ok++;
					all_diff.push(groupIdx.toString().concat(",",groupIdxRef.toString()));
				}

			}
		}
		var tot_k = groupN;
		//_log_debug(all_diff)
		//_log_debug()

		for ( var groupIdx = 0; groupIdx < groupN; groupIdx++ )
		{
			var k = groupArr[groupIdx];
			//_log_debug(k)
			var nk = 0;
			for (var j = 0; j < all_diff.length ; j++ )
			{
				var a = all_diff[j];
				if(a.indexOf(k.toString()) != -1) nk++;
			}
			if(nk == tot_k-1) {
				sig_wilc[i] = true;
				break;
			}
		}
		if(sig_wilc[i] == undefined) sig_wilc[i] = false;
	}

	var sub_biomarkers = [];
	var max_p_value = 0;
	//_log_debug(sig_wilc.length)
	for (var i = 0; i < sig_wilc.length; i++ )
	{
		if(sig_wilc[i] == true){
			sub_biomarkers.push(biomarkers[i]);
			if(max_p_value < biomarkers[i].pValue[0]) max_p_value = biomarkers[i].pValue[0];
		}
	}

	_log_debug("Number of significantly discriminative features:" + sub_biomarkers.length + "(" + biomarkers.length + ") before internal wilcoxon");
	if(sub_biomarkers.length == 0) {
		_log_debug("No features with significant differences between the two classes");
		var lefse_result = {lefse_result : [], group_names : input.groupName};
		return lefse_result;
	}


	var limit_num = 300;
	var target_num = math.min(math.max(sampleN, 100),300);
	//_log_debug("#sample: "+ sampleN)
	//_log_debug("#group: "+ groupN)
	//if(sub_biomarkers.length > target_num & sub_biomarkers.length < limit_num) _log_debug("Because the number of features (biomarkers) is over " + limit_num + ", LEfSe will select " + target_num  + " features (biomarkers) using 'p-value' to prevent the memory overflow.");
	if(sub_biomarkers.length > limit_num) _log_debug("Because the number of features (biomarkers) is over " + limit_num + ", LEfSe will select " + target_num  + " features (biomarkers) using 'p-value' to prevent the memory overflow.");

	if(sub_biomarkers.length > target_num){
		sub_biomarkers.sort(function(a,b){return a.pValue[0]-b.pValue[0];}); // sort p-value in decremental order
		//sub_biomarkers=sub_biomarkers.slice(0,limit_num-1); 
		sub_biomarkers=sub_biomarkers.slice(0,target_num);  // Because the number of features (biomarkers) is over " + limit_num + ", LEfSe will select " + target_num  + " features (biomarkers) using 'p-value' to prevent the memory overflow.
	}
	/*
	if( sub_biomarkers.length < limit_num ){
		sub_biomarkers.sort(function(a,b){return a.pValue[0]-b.pValue[0];}); // sort p-value in decremental order
		sub_biomarkers=sub_biomarkers.slice(0,target_num); // pick top_#target_number biomarkers
	}
	*/
	/*
	if(sub_biomarkers.length > target_num){
		sub_biomarkers.sort(function(a,b){return b-a;}); // sort p-value in decremental order
		sub_biomarkers=sub_biomarkers.slice(0,target_num-1); // pick top_#target_number biomarkers
	}
	*/
	/*
	while(sub_biomarkers.length > limit_num | sub_biomarkers.length > target_num & max_p_value != 0){
		var tmp_biomarkers = [];
		var second_max_p_value = 0;
		for (var i = 0; i < sub_biomarkers.length; i++)
		{
			if(sub_biomarkers[i].pValue[0] < max_p_value) tmp_biomarkers.push(sub_biomarkers[i]);
			if(max_p_value > sub_biomarkers[i].pValue[0] & second_max_p_value < sub_biomarkers[i].pValue[0]) second_max_p_value = sub_biomarkers[i].pValue[0];
		}
		sub_biomarkers = tmp_biomarkers;
		max_p_value = second_max_p_value;
	}
	*/
	//console.timeEnd("willcoxon")
	_log_debug("Now, LEfSe will execute LDA with " + sub_biomarkers.length + " features.");


// 3. Signed LDA log-score
	//console.time("matX")
	var matX = new Array(groupN);

	for (var groupIdx = 0; groupIdx < groupN; groupIdx++)
	{
		matX[groupIdx] = new Array(sub_biomarkers.length);
	}
	
	for (var i = 0; i < sub_biomarkers.length; i++)
	{
		var arr = sub_biomarkers[i].ratio.map(function(value){return value*10000;});
		var tmp_arr = [];
		for (var groupIdx = 0; groupIdx < groupN ; groupIdx++)
		{
			tmp_arr[groupIdx] = [];
		}
		for (var j = 0; j < groupInput.length ; j++){
			tmp_arr[groupInput[j]].push(arr[j]);
		}
		for (var groupIdx = 0; groupIdx < groupN ; groupIdx++)
		{
			matX[groupIdx][i] = tmp_arr[groupIdx];
		}
	}

	// override matX[groupIdx], 
	// it seems to be adjust 0 0 0 0 --> some values... prevent matrix be singular
	// using  normal(0.0,Math.max(value*0.05,0.01) makes big value much bigger....
	for (var i = 0; i < sub_biomarkers.length; i++)
	{
		for ( var groupIdx = 0; groupIdx < groupN; groupIdx++ )
		{
			var v = matX[groupIdx][i];
			if(set(v).length > Math.max(v.length * 0.5, 4)) continue;
			var v2 = v.map(function(value){
				//return math.abs(normal(0.0,Math.max(value*0.05,0.01)) + value);
				return math.abs(normal(0.0,Math.max(value*0.000000005,0.000000001))+value);
			});
			matX[groupIdx][i] = v2;
		}
	}
	
	var matX_t = [];
	var mean_arr_ld = [];
	for (var groupIdx = 0; groupIdx < groupN; groupIdx++)
	{
		mean_arr_ld.push(math.mean(matX[groupIdx],1));
		matX_t.push(math.transpose(matX[groupIdx]));
	}

	
	var fract_sample = 1; //in jm's code 0.667
	var n_boots = 1; //in lefse original program, nboots=30
	var rfk = parseInt((sampleN) * fract_sample);
	var ncl = groupN;
	var min_cl = parseInt(math.min(nByGroup) * fract_sample * fract_sample * 0.5);
	var min_cl = Math.max(min_cl,1);
	var boot_means = new Array(sub_biomarkers.length);

	for (var i = 0; i < boot_means.length ; i++ )
	{
		boot_means[i] = new Array(n_boots);
		for (var n_boot = 0; n_boot < n_boots ; n_boot++ )
		{
			boot_means[i][n_boot] = [];
		}
	}


	//*****************
	// not bootstrap

	//console.time("lda_AxB");
	for (var groupIdxA = 0; groupIdxA < (groupN-1); groupIdxA++)
	{
		for (var groupIdxB = (groupIdxA + 1); groupIdxB < groupN ; groupIdxB++)
		{
			var paired_matX = [matX[groupIdxA], matX[groupIdxB]];
			var paired_matX_t = [matX_t[groupIdxA], matX_t[groupIdxB]];


			var paired_mean_arr = [mean_arr_ld[groupIdxA], mean_arr_ld[groupIdxB]];
			var paired_sampleN = matX_t[groupIdxA].length + matX_t[groupIdxB].length;
			var prop = math.divide([matX_t[groupIdxA].length, matX_t[groupIdxB].length], paired_sampleN);

			var lda_result = lda(paired_matX,paired_matX_t,paired_mean_arr,2,paired_sampleN,sub_biomarkers,prop);
			//console.log(lda_result)
			for (var i = 0; i < sub_biomarkers.length ; i++)
			{
				boot_means[i][0].push(lda_result[i]);
			}
		}
	}

	//console.timeEnd("lda_AxB");
	
	//bootstapping
/*
	for (var n_boot = 0; n_boot < n_boots ; n_boot++ )
	{

		for (var rtmp = 0; rtmp < 1000; rtmp++)
		{
			var rand_s = [];
			for (var fk = 0; fk < rfk ; fk++)
			{
				rand_s.push(getRandomArbitrary(0, sampleN));
			}
			if(!contant_within_classes_or_few_per_class(matX_t,rand_s,min_cl,ncl,groupInput,nByGroup)) break;	
		}

		var sub_matX_t = new Array(matX_t.length);

		for (var groupIdx = 0; groupIdx < matX_t.length; groupIdx++)
		{
			sub_matX_t[groupIdx] = new Array();
		}

		for (var rand = 0; rand < rand_s.length; rand++)
		{
			var idx = rand_s[rand];
			var groupIdx = groupInput[idx];

			var s_idx = 0;
			for (var j = 0; j < idx; j++)
			{
				if(groupInput[j] == groupIdx) s_idx++;
			}
			sub_matX_t[groupIdx].push(matX_t[groupIdx][s_idx]);
		}

		var sub_mean_arr = [];
		var sub_matX = [];

		for (var groupIdx = 0; groupIdx < groupN; groupIdx++)
		{
			sub_matX.push(math.transpose(sub_matX_t[groupIdx]));
			sub_mean_arr.push(math.mean(sub_matX[groupIdx],1));
		}

		for (var groupIdxA = 0; groupIdxA < (groupN-1); groupIdxA++)
		{
			for (var groupIdxB = (groupIdxA + 1); groupIdxB < groupN ; groupIdxB++)
			{
				var paired_matX = [sub_matX[groupIdxA], sub_matX[groupIdxB]];
				//_log_debug("paired_matX: " + "["+ math.size(paired_matX) + "]");
				
				var paired_matX_t = [sub_matX_t[groupIdxA], sub_matX_t[groupIdxB]];
				//_log_debug("paired_matX_t: " + math.size(paired_matX_t) + "]");
				
				var paired_mean_arr = [sub_mean_arr[groupIdxA], sub_mean_arr[groupIdxB]];
				var paired_sampleN = sub_matX_t[groupIdxA].length + sub_matX_t[groupIdxB].length;
				var prop = math.divide([sub_matX_t[groupIdxA].length, sub_matX_t[groupIdxB].length], paired_sampleN);
				
				//console.time("lda_single");
				var lda_result = lda(paired_matX,paired_matX_t,paired_mean_arr,2,paired_sampleN,sub_biomarkers,prop);
				//console.timeEnd("lda_single");
				
				for (var i = 0; i < sub_biomarkers.length ; i++)
				{
					boot_means[i][n_boot].push(lda_result[i]);
				}
			}
		}

	}

	*/

	var all_pair_len = (groupN) * (groupN-1) / 2;

	var res = [];

	var th = 0;
	//console.log(sub_biomarkers) //same as kruskal_wallis format
	//console.log(boot_means.length)
	for (var k = 0; k < sub_biomarkers.length ; k++)
	{
		var tmp_mat = boot_means[k];
		var tmp_mat_t = math.transpose(tmp_mat);
		var tmp_arr = [];

		for (var p = 0; p < all_pair_len ; p++)
		{
			try{
				tmp_arr[p] = math.mean(tmp_mat_t[p]); // as "e" appears.... [e][r][r][o][r] --> when matrix is singular
				//_log_debug(tmp_arr[p])
			}catch(err){
				//_log_debug("matrix is singular")
				break;
			}
		}

		try{
			var m = math.max(tmp_arr);
			var eff = math.sign(m) * math.log10(math.abs(m) + 1);
			if(eff > 2.0) {
				th++;
			}
			sub_biomarkers[k]['effect_size'] = eff;
			res.push(sub_biomarkers[k]);
		}catch(err){
			
		}

	}

	_log_debug("Number of discriminative features with abs LDA score > 2.0 : " + th);
	//_log_debug(res.length)
	if(res.length == 0) {
		_log_debug("input matrix is singluar...");
		var lefse_result = {lefse_result : [], group_names : input.groupName};
		return lefse_result;
	}


	var lefse_result = {lefse_result : res, group_names : input.groupName};
	lefse_result.cladogram = _make_cladogram_data(input.samples,lefse_result);


	return lefse_result;

}

// functions

function get_within_class_variance(matX, mean_arr, sampleN){

	var diff_mat = [];

	for (var r = 0; r < matX.length; r++ )
	{
		diff_mat[r] = numeric.sub(matX[r],mean_arr[r]);
	}

	var multiplied_mat = numeric.dot(diff_mat,numeric.transpose(diff_mat));
	var var_mat = numeric.div(multiplied_mat,sampleN);

	return var_mat;

}


function get_matrix_subtract_multiply(matX, second_mat, mean_arr){
	
	var diff_mat = [];

	for (var r = 0; r < matX.length ; r++ )
	{
		diff_mat[r] = numeric.sub(matX[r],mean_arr[r]);
	}

	var multiplied_mat = numeric.dot(numeric.transpose(diff_mat), second_mat);

	return multiplied_mat;

}

function set (arr) {
	return arr.reduce(function (a, val) {
		if (a.indexOf(val) === -1) {
			a.push(val);
		}
		return a;
	}, []);
}

function lda(matX,matX_t,mean_arr,groupN,sampleN,sub_biomarkers,prop){


	var tol = 0.0000000001; // A tolerance to decide if a matrix is singular
	var within_v = [];

	for (var groupIdx = 0; groupIdx < groupN ; groupIdx++)
	{
		within_v[groupIdx] = get_within_class_variance(matX[groupIdx], mean_arr[groupIdx], sampleN-1);
	}
	var sum_within_v = within_v[0];

	for (var groupIdx = 1; groupIdx < groupN ; groupIdx++ )
	{
		sum_within_v = numeric.add(sum_within_v, within_v[groupIdx]);
	}

	var diag_mat = math.diag(sum_within_v);
	var f1 = math.sqrt(diag_mat);

	// checking variables to be constant within groups
	for (var i = 0; i < f1.length ; i++)
	{
		if(f1[i] < tol){
			_log_debug("CANNOT OPERATE LDA (singular matrix)");
			return "error";
		}
		f1[i] = 1/f1[i];
	}

	// scale columns to unit variance before checking for collinearity : mle method
	var scaling = math.diag(f1);

	var fac = 1 / sampleN;
	
	var X = [];
	var tmpX;

	for (var groupIdx = 0; groupIdx < groupN; groupIdx++) // merged Xs
	{
		tmpX = get_matrix_subtract_multiply(matX[groupIdx], scaling, mean_arr[groupIdx]);

		for (var i = 0; i < tmpX.length ; i++ )
		{
			X.push(math.multiply(tmpX[i],math.sqrt(fac)));
		}		
	}
	
	// Need more rows than columns, So do not transpose later, use V in numeric library
	if(math.transpose(X).length < X.length){
		var X_s = numeric.svd(X);
	} else {
		var X_s = numeric.svd(math.transpose(X));
	}

	var rank = 0;
	for (var i = 0; i < X_s.S.length ; i++ )
	{
		if(X_s.S[i] > tol) rank++;
	}
	
	if(rank == 0) return "rank = 0: variables are numerically constant";
	if(rank < sampleN) _log_debug("variables are collinear");

	var invS = [];
	var newV = [];
	for (var i = 0; i < sub_biomarkers.length ; i++)
	{
		var tmp_arr = [];
		for (var j = 0; j < rank; j++)
		{
			if(numeric.transpose(X).length < X.length){
				tmp_arr.push(X_s.V[i][j]);
			} else {
				tmp_arr.push(X_s.U[i][j]);
			}
		}
		newV.push(tmp_arr);
	}
	for (var j = 0; j < rank; j++)
	{
			invS.push(1/X_s.S[j]);
	}

	var diag_mat_S = math.diag(invS);
	var new_scaling = numeric.dot(newV,diag_mat_S);
	new_scaling = math.multiply(scaling,new_scaling);

	var xvar = math.multiply(math.transpose(prop), mean_arr);
	
	var fac = 1/groupN;

	var tmpX = [];
	for (var groupIdx = 0; groupIdx < groupN; groupIdx++)
	{
		tmpX[groupIdx] = numeric.sub(mean_arr[groupIdx],xvar);
	}
	//_log_debug(math.size(tmpX));
	//_log_debug(new_scaling);
	tmpX = math.multiply(tmpX,new_scaling);
	var fac = math.sqrt(math.multiply(math.multiply(prop,sampleN),fac));
	
	var X = [];
	for (var groupIdx = 0; groupIdx < groupN; groupIdx++)
	{
		X[groupIdx] = math.multiply(tmpX[groupIdx], fac[groupIdx]);
	}

	var X_s = numeric.svd(numeric.transpose(X));
	
	var rank = 0;
	for (var i = 0; i < X_s.S.length ; i++ )
	{
		if(X_s.S[i] > tol * X_s.S[0]) rank++;
	}

	if(rank == 0) return "rank = 0: variables are numerically constant";

	var X_s_tU = numeric.transpose(X_s.U);

	var final_scaling = [];
	for (var i = 0; i < sub_biomarkers.length ; i++)
	{
		final_scaling.push(math.multiply(new_scaling[i],X_s_tU[0]));
	}

	var w_unit = math.divide(final_scaling,math.sqrt(math.sum(math.dotPow(final_scaling,2))));
	//var LD = [];
	var effect_sizes =[];
	for (var groupIdx = 0; groupIdx < groupN; groupIdx++)
	{
		var tmp_arr = [];
		for (var i = 0; i < matX_t[groupIdx].length ; i++)
		{
			//LD.push(numeric.dot(matX_t[groupIdx][i],w_unit));
			tmp_arr.push(numeric.dot(matX_t[groupIdx][i],w_unit));
		}
		effect_sizes[groupIdx] = math.mean(tmp_arr);
	}
	var effect_size = math.abs(numeric.sub(effect_sizes[0],effect_sizes[1]));
	var wfinal = math.abs(numeric.dot(w_unit,effect_size));

	var gm = math.abs(numeric.sub(mean_arr[0], mean_arr[1]));
	var res = [];

	for (var i = 0; i < gm.length; i++)
	{
		res.push((gm[i] + wfinal[i]) * 0.5);
	}
	
	return res;

}

function permanova(input,server_log_debug){
		if (server_log_debug != undefined && server_log_debug != null) {
			_log_debug = server_log_debug;
		}

		if(input.group.unique().length < 2){
			_log_debug("#group must more than 2 in permanova");
		}else{
			_log_debug("[info] operate permanova operation for "+input.distance_matrix.length+" distance metrics...")
			let result = []

			//for(let n_dist = 0; n_dist<1;n_dist++){ / for test
			for(let n_dist = 0; n_dist<input.distance_matrix.length;n_dist++){
				_log_debug("[info] calculate p-value for all("+input.group.unique().length+") groups")
				let main_stat = [];
				let pair_upper = []; 
				let group = input.group
				let matrix = input.distance_matrix[n_dist]
				let d2 = numeric.pow(matrix,2)
							
				let P_tot = permanova_core(group.slice(),d2,999);			
				main_stat.push(P_tot)

				_log_debug("[info] calculate inner-statistics...")
				let inner = []
				for(let i=0;i<input.group.unique().length;i++){
					let inner_group_check  = group.slice().map(function(v){if(v==i)return 1; else return 0;})
					let inner_group = numeric.mul(inner_group_check,input.group.slice().map(v =>v+1))
					
					let adjust = numeric.mul(numeric.identity(inner_group.length),0.00000000000001)
					let inner_d2 = numeric.mul(numeric.dot(numeric.transpose([inner_group_check]),[inner_group_check]),numeric.add(d2,adjust))
					
					let inner_group_red = inner_group.filter(elem => elem !==0).map(v => v-1)
					let adjust_red = numeric.mul(numeric.identity(inner_group_red.length),0.00000000000001)
					let inner_d2_red = numeric.sub(inner_d2.map(v => v.filter(elem => elem !==0)).filter(elem => elem.length !== 0),adjust_red)
					
					let inner_stat = {F:null,p:null,key:0}
					let basic = stat_for_boxplot(inner_d2_red);
					for(let attrname in basic){inner_stat[attrname] = basic[attrname]}

					inner.push(inner_stat)

				}

				if(input.group.unique().length > 1){
					_log_debug("[info] calculate pairwise p-value for "+0.5*input.group.unique().length*(input.group.unique().length-1) + " permutations")
					//split pairwise_group
					let upper_idx = 0
					let pair_tmp = []
					for(let i=0;i<input.group.unique().length-1;i++){
						for(let j=i+1;j<input.group.unique().length;j++){
							//prune pair_group
							let pair_group_check = group.slice().map(function(v){if(v==i || v==j) return 1; else return 0})
							let pair_group = numeric.mul(pair_group_check,input.group.slice().map(v => v+1))
							//console.log(pair_group.length)
							//let adjust = numeric.mul(numeric.identity(pair_group.length),0.000000000000000000001)
							let adjust = math.multiply(0.00000000000001,math.ones(pair_group.length,pair_group.length))
							//console.log(adjust)
							let pair_d2 = numeric.mul(numeric.dot(numeric.transpose([pair_group_check]),[pair_group_check]),numeric.add(d2,adjust._data)) //couldn't reduce matrix size....for sake of time....
							//console.log(pair_d2)
							//reduce matrix size
							let pair_group_red = pair_group.filter(elem => elem !==0).map(v => v-1)
							//let adjust_red = numeric.mul(numeric.identity(pair_group_red.length),0.000000000000000000001)
							let adjust_red = math.multiply(0.00000000000001,math.ones(pair_group_red.length,pair_group_red.length))
							let pair_d2_red = numeric.sub(pair_d2.map(v => v.filter(elem => elem !==0)).filter(elem => elem.length 	!== 0),adjust_red._data)
							//console.log(adjust_red)

							let stat = permanova_core(pair_group_red,pair_d2_red,999)
							stat.key = upper_idx
							
							upper_idx++
							let basic = stat_for_boxplot(pair_d2_red);
							for(let attrname in basic){stat[attrname] = basic[attrname]}

							//console.log(stat)
							pair_tmp.push(stat)
							//so pair_stat holds upper part of p-value matrix...
						}
					}


					//calculate fdr
					getQ(pair_tmp);
					//make pair stat into matrix

					let pair_idx= 0;
					let pair_stat=math.zeros(input.group.unique().length,input.group.unique().length)._data
					for(let i=0;i<input.group.unique().length-1;i++){
						for(let j=i+1;j<input.group.unique().length;j++){
							pair_stat[i][j] = pair_tmp[pair_idx]
							pair_stat[j][i] = pair_tmp[pair_idx]
							pair_idx++
						}
					}
					for(let i=0;i<input.group.unique().length;i++){
						pair_stat[i][i] = inner[i];
					}
					pair_upper.push(pair_stat)
					//console.log(pair_stat)


				}else{
					_log_debug("[info] #group less than 2, return result!")
					//P_tot.q = P_tot.p
					//console.log(P_tot)
					//let pair_idx= 0;
					//let pair_stat=math.zeros(input.group.unique().length,input.group.unique().length)._data //2x2
					//pair_stat[0][1]=P_tot
					//pair_stat[1][0]=P_tot
					//for(let i=0;i<input.group.unique().length;i++){
					//	pair_stat[i][i] = inner[i];
					//}

					pair_upper.push([])
					
				}

				//let res = {main_stat: main_stat, pair_stat:pair_upper,distance_matrix:input.distance_matrix[n_dist]}
				let res = {main_stat:main_stat,pair_stat: pair_upper}

				result.push(res)


			}
		let Result = {permanova_result: result, group: input.group, groupName: input.groupName}
		return Result
	}
		
}


function picrust(input){
        let group = input.group
        let profile = input.functional_profiles
        let pV = []
        for(let i=0; i<profile.length; i++){
                //split by group
                let in_split = profile[i].v.chunk(wordcount(group).count)
                pV.push({name : profile[i].nm, pvalue : kruskal_core.apply(this,in_split).pvalue, profile : in_split});
        }
        let pV_filter = pV.filter(function(v){return v.pvalue<=0.05})
        //console.log(pV_filter[1].profile)
        let matX = []
        for(let g=0;g<group.unique().length;g++){
        	let tmp = []
	        for(let i=0;i<pV_filter.length;i++){
        		tmp.push(pV_filter[i].profile[g])
        		//matX[1][i] = pV_filter[i].profile[1]
        	}
        	matX.push(tmp)
        }
        //console.log(math.size(matX[0]))
        for(let i=0; i<group.unique().length;i++){
        	for(let j=0;j<pV_filter.length;j++){
        		//console.log(i+" "+j)
        		matX[i][j] = matX[i][j].map(function(value){return math.abs(normal(0.0,Math.max(value*0.000000005,0.000000001))+value)}) //to prevent singular...?
        	}
        }
        //console.log(math.size(matX))
        //console.log(matX)
        let lda_result = []
        for(let gIdx1=0;gIdx1<group.unique().length-1;gIdx1++){
        	for(let gIdx2=gIdx1=1; gIdx2<group.unique().length;gIdx2++){
        		let lda_result_tmp = lda_new(matX[gIdx1],matX[gIdx2])
        		lda_result.push(lda_result_tmp)
        		//console.log(lda_result)
        	}
        }


        for(let k=0;k<pV_filter.length; k++){ //when it comes to multiple-group, call//
        	let m = lda_result[k];
        	pV[k].effect_size = math.sign(m) * math.log10(math.abs(m)+1)
        }
        return pV


}


function lda_new(input){
	//matX : m x d, matY : n x d
	const args =  Array.from(arguments);
	if(args.length !== 2){
		console.error("[E] lda need 2 input");
	}else{
		const arr = args.reduce(function(p,c){return p.concat(c)});

		const tol = 0.0000000001
		let arr_split = arr.slice().chunk(args.map(function(v){return v.length}));

		let m1 = math.mean(arr_split[0],1)
		console.log("???????????????????????????")
		console.log(math.size(arr_split[0]))
		console.log(math.size(arr_split[1]))
		console.log("???????????????????????????")
		let m2 = math.mean(arr_split[1],1)
		let m = math.mean(arr,1)

		let diff1 = arr_split[0].map(function(v,i){return numeric.sub(v,m1[i])});
		let diff2 = arr_split[1].map(function(v,i){return numeric.sub(v,m2[i])});


		let _checker = numeric.svd(numeric.transpose(numeric.transpose(diff1).concat(numeric.transpose(diff2))))
		/*
		let new_diff ; 
		if(_checker.S.some(function(elem){return elem <= tol})){
			console.log("!!")
			let rank = _checker.S.filter(function(elem){return elem >tol}).length;
			if(rank == 0 ){
				console.error("rank = 0: variables are numerically constant ")
			}else{
				console.info("varialbes are collinear")
				let newU = _checker.U.map(function(v){return v.slice(0,rank)});
				let newS = _checker.S.slice(0,rank);
				let newV = _checker.V.map(function(v){return v.slice(0,rank)}).slice(0,rank);
				

				console.log("U: "+math.size(newU))
				console.log("S: "+math.size(newS))
				console.log("V: "+math.size(newV))
				

				//S_W_new = numeric.dot(numeric.dot(newU,newS),numeric.transpose(newV));
				//console.log()
		
				new_checker = numeric.dot(numeric.dot(newU,[newS]),newV);

			}
		}else{
			console.log("same rank")
			new_diff = numeric.transpose(numeric.transpose(diff1).concat(numeric.transpose(diff2)));
		}

		let new_diff_split = numeric.transpose(new_diff).slice().chunk(args.map(function(v){return v.length}));
		*/



		let sigma_1 = numeric.div(numeric.dot(diff1,numeric.transpose(diff1)),arr_split[0][0].length-1)
		let sigma_2 =numeric.div(numeric.dot(diff2,numeric.transpose(diff2)),arr_split[1][0].length-1)
		let S_B = numeric.dot(numeric.transpose([numeric.sub(m1,m2)]),[numeric.sub(m1,m2)]);
		let S_W = numeric.add(sigma_1,sigma_2);
		console.log("TTTTTTTTTTTTTT")
		console.log(math.size(S_B))
		console.log("TTTTTTTTTTTTTT")


		console.log(math.size(diff1))
		console.log(math.size(diff2))
		//console.log(math.size(numeric.transpose(diff1).concat(numeric.transpose(diff2))))
		//let _checker = numeric.svd(numeric.transpose(numeric.transpose(diff1).concat(numeric.transpose(diff2))))
		console.log("xxxxxxxxxxxxxxxx")
		console.log(math.size(_checker.U))
		console.log(math.size(_checker.S))
		console.log(math.size(_checker.V))
		console.log("xxxxxxxxxxxxxxxx")
		//console.log(_checker.S)

				/*
		console.log(math.size(S_W))
		let S_W_new;
		let _svd = numeric.svd(S_W);
		let _eig = numeric.eig(S_W)
		//console.log(_eig)
		if(_svd.S.some(function(elem){return elem <= tol})){
			console.log("!!")
			let rank = _svd.S.filter(function(elem){return elem >tol}).length;
			if(rank == 0 ){
				console.error("rank = 0: variables are numerically constant ")
			}else{
				console.info("varialbes are collinear")
				let newS = _svd.S.slice(0,rank);
				let newV = _svd.V.slice(0,rank);
				let newU = _svd.U.map(function(v){return v.slice(0,rank)});

				console.log("U: "+math.size(newU))
				console.log("S: "+math.size(newS))
				console.log("V: "+math.size(newV))
				

				//S_W_new = numeric.dot(numeric.dot(newU,newS),numeric.transpose(newV));
				//console.log()
		
				S_W_new = numeric.dot(numeric.dot(newU,[newS]),newV);

			}
		}else{
			S_W_new = S_W;
		}
		*/
		let S_W_new = S_W;
		console.log(numeric.inv(S_W_new))
		console.log(numeric.svd(numeric.dot(numeric.inv(S_W_new),S_B)))
		let w_star = numeric.abs(numeric.svd(numeric.dot(numeric.inv(S_W_new),S_B)).U.map(function(v){return v[0]}));
		console.log(w_star)
		console.log(math.size(w_star))
		let x0_transform = numeric.dot(numeric.transpose(arr_split[0]),w_star)
		console.log(math.size(x0_transform))
		console.log(x0_transform)
		let x1_transform = numeric.dot(numeric.transpose(arr_split[1]),w_star)
		console.log(math.size(x1_transform))

		let eff = math.abs(numeric.sub(math.mean(x0_transform),math.mean(x1_transform)))
		//console.log(eff)
		let w_final = math.abs(numeric.dot(w_star,eff));

		console.log("::::::::::::::::")
		//console.log( w_final )
		//console.log(math.size(w_final))
		console.log("::::::::::::::::")
		return numeric.mul(0.5,numeric.add(math.abs(numeric.sub(m1,m2)),w_final))

	}
}


function stat_for_boxplot(matrix){
	//
	//console.log(matrix)
	let d = numeric.sqrt(matrix)
	let res ={}
	let tmp = []
	//get upper part of matrix
	for(let i=0;i<matrix.length-1;i++){
		for(let j=i+1;j<matrix.length;j++){
			tmp.push(matrix[i][j])
		}
	}
	let tmp_sort = tmp.slice().sort(function(a,b){return a-b})
	console.log(math.size(d))
	//console.log(tmp)
	res.median = math.median(tmp_sort)
	//res.q1 = math.quantileSeq(tmp,1/4,true) // mathjs interpolates quantile... 
	//res.q3 = math.quantileSeq(tmp,3/4,true) // mathjs interpolates quantile...
	res.q1 = tmp_sort[Math.floor(tmp.length/4)]
	res.q3 = tmp_sort[3*Math.floor(tmp.length/4)]
	res.max = math.max(tmp_sort)
	res.min = math.min(tmp_sort)
	res.avg = math.mean(tmp_sort)
	res.std = math.std(tmp_sort)
	res.var = math.var(tmp_sort)
	console.log(res)
	return res
}

function permanova_core(group,matrix,n_permutation){
	if(group.unique().length < 2){
			_log_debug("#group must more than 2 in permanova");
		}else{
			let F_list = [];
			let F_ori = getF(group,matrix);

			for(let i=0;i<n_permutation;i++){
				let tmp = group.slice()
				shuffle(tmp)
				F_list.push(getF(tmp,matrix))
			}

			let P = (F_list.filter(function(elem){return elem > F_ori}).length + 1)/(n_permutation+1)

			return {F : F_ori, p : P}
	}
}


function getF(group,matrix){
	if(group.unique().length < 2){
		_log_debug("#group must more than 2 to get F-value");
	}else{
		let sst = numeric.sum(matrix) * 0.5 //same speed, but shorter code
		let ssw = get_ssw(group,matrix);
		let ssa = (1/group.length)*sst-ssw	
		let F = (ssa/(wordcount(group).word.length-1)) / (ssw/(group.length-wordcount(group).word.length))
		return F
	}
}

function get_ssw(group,matrix){
	let ssw = 0
	for(let i=0;i<group.unique().length;i++){
		let checkG = group.map(function(v){if(v==group.unique()[i]) return 1; else return 0})
		ssw += 0.5*numeric.sum(numeric.mul(numeric.dot(numeric.transpose([checkG]),[checkG]),matrix)) * (1/checkG.filter(function(elem){return elem ==1}).length)
	}
	return ssw;
}

function shuffle(input){
	for(let i=input.length; i; i-=1 ){
		let j = Math.floor(Math.random()*i);
		let x = input[i-1]
		input[i-1] = input[j];
		input[j] = x
	}
}





module.exports = {
	biomarker_xor: xor,
	biomarker_wilcoxon: wilcoxon,
	biomarker_kruskalWallis: kruskal,
	biomarker_lefse: lefse,
	xor_core : xor_core,
	wilcoxon_core : wilcoxon_core,
	kruskal_core : kruskal_core,
	permanova : permanova,
	permanova_core : permanova_core,
	statistic : stat_for_boxplot,
	picrust : picrust

}