// High Luminance Color Appearance Model Extractor
// 
// This javascript program will extract all the experiemntal data from the file data_raw.txt
// The text file is a copy and past of the PDF file found at http://jankautz.com/publications/experimental_data.pdf
// 
// The extrator will create two files 
// 1. data_phase.csv - holds all the viewing conditions data for each of the phases
// 2. data_perceptual.csv - holds all the physical measurements and perceptual estimates merged with viewing conditions 


const fs = require('fs');

const text = fs.readFileSync(`data_raw.txt`, 'utf-8');
const phaseData = extractPhaseData(text);
console.log(`Extracted ${phaseData.split('\n').length} for data_phase.csv`);
fs.writeFileSync(`data_phase.csv`, phaseData);

const allData = extractAllData(text, phaseData);
console.log(`Extracted ${allData.split('\n').length} for data_perceptual.csv`);
fs.writeFileSync(`data_perceptual.csv`, allData);


function extractPhaseData(data) {
	return data
		.replace(/Table 1: Summary of viewing conditions for all 19 phases[\s\S]*/gm, '')
		.replace(/[\s\S]*Luminance/gm, 'Phase Medium CCT Xw Yw Zw La Pb Xb Yb Zb Ambient')
		.replace(/%/gm,'')
		.replace(/ /gm, ',')
	;
}

function extractAllData(data, phaseData) {
	const phaseDataSplit = phaseData.split('\n');
	let result = 'Phase Medium CCT Xw Yw Zw La Pb Xb Yb Zb Ambient Color Xs Ys Zs Js Ms Hs\n';
	for (var i = 1; i < 20; i++) {
		result += extractData(text, i, phaseDataSplit);
	}
	return result
		.replace(/^$\n/gm,'')
		.replace(/ /gm, ',')
		.replace(/N\/A/gm,'NA')
	;
}

function extractData(data, n, phaseDataSplit) {

	// get the Table no to extract (all data is paired up except for Phase 19) 
	const phaseTableNo = ((n % 2 === 1) && (n < 19))? n+1 : n; 

	// extract the table
	const table = data
		.replace(RegExp(`[\\s\\S]*Phase ${phaseTableNo}\nColor`,'gm'), '')
		.replace(/Table [\s\S]*/gm,'')
		.replace(/^[^\d].*/gm,'')
	;
	const halfTable = (n % 2 === 1) ?
		table.replace(/(^\d+ \S+ \S+ \S+ \S+ \S+ \S+).*/gm, '$1') :
		table.replace(/(^\d+) \S+ \S+ \S+ \S+ \S+ \S+ (\S+ \S+ \S+ \S+ \S+ \S+).*/gm, '$1 $2');

	return halfTable
		.replace(/(^\S+) (.*)/gm, `${phaseDataSplit[n]} $1 $2`)
		.replace(/ /gm, ',')
	;
}


// takes a string or array of strings and normalises it
function normalizeString(str) {
	return (str)
		.toLowerCase()
		// remove punctuation marks
		.replace(/[,\.!:;\?\(\)\[\]\{\}]/g, '')
		.replace(/['’]/g, '')			// fast match for common appostrophe
		.replace(/"/g, '')			    // fast match for common double quote
		.replace(/-/g, ' ')			    // fast match for common dash
		// coerce single curly quotes
		.replace(/[\u2018\u2019\u201A\u201B\u2032\u2035]+/g, '')
		// coerce double curly quotes
		.replace(/[\u201C\u201D\u201E\u201F\u2033\u2036]+/g, '');
}
