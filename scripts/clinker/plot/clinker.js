function serialise(svg) {
	/* Saves the figure to SVG in its current state.
	 * Clones the provided SVG and sets the width/height of the clone to the
	 * bounding box of the original SVG. Thus, downloaded figures will be sized
	 * correctly.
	 * This function returns a new Blob, which can then be downloaded.
	*/
	let node = svg.node()
		.cloneNode(true)

	// Remove all hidden elements
	d3.select(node)
		.selectAll(".hidden")
		.remove()

	// Draw filtered node in hidden <div> to calculate real bounding box
	let container = svg.node().parentNode
	let div = d3.select(container)
		.append("div")
		.attr("visibility", "hidden")
	div.append(() => node)
	let bbox = div.select("g.clusterMapG").node().getBBox()
	div.remove()

	const xmlns = "http://www.w3.org/2000/xmlns/"
	const xlinkns = "http://www.w3.org/1999/xlink"
	const xhtml = "http://www.w3.org/1999/xhtml"
	const svgns = "http://www.w3.org/2000/svg"

	node.setAttribute("width", bbox.width);
	node.setAttribute("height", bbox.height);
	node.setAttributeNS(xmlns, "xmlns", svgns);
	node.setAttributeNS(xmlns, "xmlns:xlink", xlinkns);
	node.setAttributeNS(xmlns, "xmlns:xhtml", xhtml);

	// Adjust x/y of <g> to account for axis/title position.
	// Replaces the transform attribute, so drag/zoom is ignored.
	d3.select(node)
		.select("g.clusterMapG")
		.attr("transform", `translate(${Math.abs(bbox.x)}, ${Math.abs(bbox.y)})`)

	const serializer = new window.XMLSerializer;
	const string = serializer.serializeToString(node);
	return new Blob([string], {type: "image/svg+xml"});
}

function gnat_serialise(svg) {
	/* Saves the figure to SVG in its current state.
   * Without the cluster names and legend, only the clusters and genes are shown.
   * Called by the GNAT app to save each cluster as a SVG for the dashboard and previews.
	 * Clones the provided SVG and sets the width/height of the clone to the
	 * bounding box of the original SVG. Thus, downloaded figures will be sized
	 * correctly.
	 * This function returns a new Blob, which can then be downloaded.
	*/

	let node = svg.node()
		.cloneNode(true)

	// Remove all hidden elements
	d3.select(node)
		.selectAll(".hidden")
		.remove()

  // Remove all <g> elements in the legend
  d3.select(node)
    .selectAll("g.legend")
    .remove()

  // Remove all <g> elements in the scale bar
  d3.select(node)
    .selectAll("g.scaleBar")
    .remove()

  // Remove all <g> elements in the colour bar
  d3.select(node)
    .selectAll("g.colourBar")
    .remove()

  // Remove all <g> elements in clusters with class clusterInfo
  d3.select(node)
    .selectAll("g.clusterInfo")
    .remove()

	// Draw filtered node in hidden <div> to calculate real bounding box
	let container = svg.node().parentNode
	let div = d3.select(container)
		.append("div")
		.attr("visibility", "hidden")
	div.append(() => node)
	let bbox = div.select("g.clusterMapG").node().getBBox()
	div.remove()

	const xmlns = "http://www.w3.org/2000/xmlns/"
	const xlinkns = "http://www.w3.org/1999/xlink"
	const xhtml = "http://www.w3.org/1999/xhtml"
	const svgns = "http://www.w3.org/2000/svg"

	node.setAttribute("width", bbox.width);
	node.setAttribute("height", bbox.height);
	node.setAttributeNS(xmlns, "xmlns", svgns);
	node.setAttributeNS(xmlns, "xmlns:xlink", xlinkns);
	node.setAttributeNS(xmlns, "xmlns:xhtml", xhtml);

	// Adjust x/y of <g> to account for axis/title position.
	// Replaces the transform attribute, so drag/zoom is ignored.
	d3.select(node)
		.select("g.clusterMapG")
		.attr("transform", `translate(${Math.abs(bbox.x)}, ${Math.abs(bbox.y)})`)

	const serializer = new window.XMLSerializer;
	const string = serializer.serializeToString(node);
	return new Blob([string], {type: "image/svg+xml"});
}

function download(blob, filename) {
	/* Downloads a given blob to filename.
	 * This function appends a new anchor to the document, which points to the
	 * supplied blob. The anchor.click() method is called to trigger the download,
	 * then the anchor is removed.
	*/
	const link = document.createElement("a");
	link.href = URL.createObjectURL(blob);
	link.download = filename;
	document.body.appendChild(link);
	link.click();
	document.body.removeChild(link);
}

//-----------------------------------------------------------------------------------------------------------------------
// This block is a customization added to the original Clinker code

let flipLoci = [];
let matchGenes = {};
let geneUidMap = {};
const numClusters = data.clusters.length;
var anchorGenes = []

for (let c in data.clusters) {
    const cluster = data.clusters[c];
    const clusterName = cluster.name.split("_");
    const eValue = clusterName[clusterName.length - 1]; // E value of psi-blast match
    const matchAccession = clusterName[clusterName.length - 2]; // Accession number of psi-blast match
    let newName = "";

    for (let n in clusterName) {
      if (n < clusterName.length - 3) {
        newName += clusterName[n] + " ";
      } else if (n == clusterName.length - 3) {
        if (clusterName[n] == 1) {
          newName += "(" + clusterName[n] + " duplicate)"
        } else {
          newName += "(" + clusterName[n] + " duplicates)"
        }
      }
    }
    data.clusters[c].name = newName;


    (cluster.loci.length > 1) ? alert(clusterName + " " + cluster.loci.length) : {}; // Alert if more than one locus per cluster.. things may not work as expected
    // Even though there should always be 1 locus passed in by the pipeline, iterate through loci[] just in case
    let flipFlag = false;
    for (let l in cluster.loci) {
        for (let g in cluster.loci[l].genes) {
            //if (cluster.name.includes("IFB5605")) {alert(cluster.loci[l].genes[g].label)}
            const gene = cluster.loci[l].genes[g];

            geneUidMap[gene.uid] = [gene.names.product, gene.label];

            if (gene.names.psi_blast_match != null) {      // If the gene is a psiblast match, append it to a dictionary to be used later for group labeling
              const matchGene = gene.names.psi_blast_match.split("(")[1].split(" ")[0];
              //matchGenes[matchGene] == null ? matchGenes[matchGene] = [gene.uid] : matchGenes[matchGene].push(gene.uid);
              matchGenes[gene.uid] = matchGene;
            }

            if (gene.label == matchAccession) {
                anchorGenes.push("gene_"+gene.uid);
                //(l == 0 && c == 0) ? anchorGeneId="gene_"+gene.uid : {} // set anchor gene if first cluster

                if (gene.strand == -1) { // Match gene is reversed, flip locus
                    flipFlag = true;
                    //break;
                }
            }
        }
        if (flipFlag) {
            // Add locus id to array to flip later
            flipLoci.push("locus_"+cluster.loci[l].uid);
        }
    }
}

// Set the additional group attributes we need for interaction, then display the 2 products
setGroupAttributes();
updateGroupProductsDisplay(document.getElementById("input-group-label-length").value);

//-----------------------------------------------------------------------------------------------------------------------


function plot(data) {
  const div = d3.select("#plot")
  const chart = ClusterMap.ClusterMap()
    .config({
      scaleFactor: 50,  // Modified, used to be 30
      cluster: {
        spacing: 50,
        alignLabels: true,
      },
      gene: {
        label: {
          show: false,
        }
      },
    })

  let plot = d3.select("#plot")
		.data([data])
    .call(chart)

  d3.select("#btn-save-svg")
    .on("click", () => {
			const svg = div.select("svg.clusterMap")
			const blob = serialise(svg)
			download(blob, "clinker.svg")
		})

	/**
	 * Saves current chart config and data to JSON.
	 * Can be loaded back into the clinker webapp using loadData().
	 */
	function saveData() {
		const obj = { config: chart.config(), data: chart.data() }
		const blob = new Blob(
			[JSON.stringify(obj, null, 2)],
			{type: "application/json"}
		);
		download(blob, "clinker.json")
	}

	/**
	 * Loads JSON generated by saveData() into clinker webapp.
	 * Should contain config and data properties.
	 */
	function loadData(event) {
		const files = event.target.files
		const file = files[0]
		if (!file) return
		const reader = new FileReader()
		reader.readAsText(file)
		reader.onload = function() {
			let newData = JSON.parse(reader.result)
			chart.config(newData.config)
			plot.selectAll("*").remove()
			plot.data([newData.data])
			plot.call(chart)
		}
		reader.onerror = function() {
			console.log("Failed to load data", reader.error)
		}
	}

	// Bind save/load
	d3.select("#btn-save-data")
		.on("click", () => saveData("clinker.json"))
	d3.select("#input-load-data")
		.on("change", loadData)

  function update(opts) {
    chart.config(opts)
    plot.call(chart)
  }

  // Populate new class label opts
  let classNames = new Set();
  if (data.classes != null) {
    if (data.classes.length > 0) {
      for (let c of data.classes) {
        // Verify the new class doesn't contain all non-labeled groups (e.g., verify every group name isn't "Group #")
        if (c.groups.length > 0) {
          let allNonLabeled = true;
          for (let g of c.groups) {
            if (!g.label.startsWith("Group ") && !(g.colour == "#808080")) {
              allNonLabeled = false;
              break;
            }
          }

          if (!allNonLabeled) {
            classNames.add(c.name);
          }
        }
      }
    }
  }

  for (let name of classNames) {
    if (name !== "Default") {
      d3.select('#select-label-class').append("option").text(name);
    }
  }

  // Add event listener to change class labels
  d3.select("#select-label-class").on("change", function() {
    let e = document.getElementById("select-label-class");
    let type = e.options[e.selectedIndex].value;

    // Check data for the class label
    let groupIdx = -1;
    for (let i = 0; i < data.classes.length; i++) {
      if (data.classes[i].name == type) {
        groupIdx = i;
        break;
      }
    }

    if (groupIdx == -1) {
      alert("Error: Class label not found in data. Please ensure the class label is present in the data.");
      return;
    } else {

      resetBolding(); // Reset bolding of groups if they were toggled on

      for (let i = 0; i < data.groups.length; i++) {
        data.groups[i].hidden = true;
      }

      // Wipe data.groups
      data.groups = [];

      // // Create a set of the added group labels
      // let addedLabeledGroups = new Set();
      // let addedLabeledGroupsColours = {};

      // Add each group from the selected class (index) to data.groups
      for (let i = 0; i < data.classes[groupIdx].groups.length; i++) {
        data.classes[groupIdx].groups[i].hidden = false;
        data.groups.push(data.classes[groupIdx].groups[i]);
      }

      // Set group nums and products
      setGroupNums();
      setGroupProducts();
    }

    // Get the number of products being displayed
    let numProducts = document.getElementById("input-group-label-length").value;
    updateGroupProductsDisplay(numProducts);

    plot.call(chart);

    // Now that the groups are added and coloured if their colours were null, go through them again to make sure groups with the same label are the same colour
    let labeledGroups = new Set();
    let labeledGroupsColours = {};
    let updatedGroups = [];

    for (let i = 0; i < data.groups.length; i++) {
      let g = data.groups[i];

      if (g.labeled) {
        // Get the core of the label (e.g., "AbiD (Homologous Group 1)" -> "AbiD")
        let coreLabel = g.baseLabel.split(" (")[0].trim();

        // Check if the group label is already in the set
        if (labeledGroups.has(coreLabel)) {
          // If it is, set this group's colour to the same as the existing one
          // alert(coreLabel + " already in set.");
          // alert("Old colour: " + g.colour + "\nNew colour: " + labeledGroupsColours[coreLabel]);
          g.colour = labeledGroupsColours[coreLabel];
        } else {
          // If it isn't, add the group label to the set and set the colour
          // alert(coreLabel + " not in set.");
          labeledGroups.add(coreLabel);
          labeledGroupsColours[coreLabel] = g.colour;
        }
      }
      
      updatedGroups.push(g);
    }

    function tokenize(label) {
      // Split label on whitespace and '_'
      return label.split(/[\s_]+/);
    }

    function isSublist(shorter, longer) {
      for (let i = 0; i <= longer.length - shorter.length; i++) {
        let match = true;
        for (let j = 0; j < shorter.length; j++) {
          if (longer[i + j] !== shorter[j]) {
            match = false;
            break;
          }
        }
        if (match) return true;
      }
      return false;
    }

    function longestCommonRun(list1, list2) {
      let maxRun = 0;

      for (let i = 0; i < list1.length; i++) {
        for (let j = 0; j < list2.length; j++) {
          let runLength = 0;
          while (
            i + runLength < list1.length &&
            j + runLength < list2.length &&
            list1[i + runLength] === list2[j + runLength]
          ) {
            runLength++;
          }
          if (runLength > maxRun) {
            maxRun = runLength;
          }
        }
      }

      return maxRun;
    }

    // Go through the groups again and parse the PSI-BLAST match groups if they haven't been coloured yet
    for (let i = 0; i < updatedGroups.length; i++) {
      let g = updatedGroups[i];
      // If the group is a PSI-BLAST match group, check if it's already been re-coloured
      if (g.match == true) {
        if (g.colour == "#808080") {
          // If it hasn't been coloured, get the PSI-BLAST match label from the top of products
          let matchLabel = g.products[0];
          let matchLabelList = tokenize(matchLabel);

          // Check if the match label is in labeledGroups (either == or as a substring or a superstring)
          let matches = [];

          for (let label of labeledGroups) {
            console.log(label)
            console.log(matchLabel)
            let labelList = tokenize(label);
            let isMatch = false;
            if (labelList.length < matchLabelList.length) {
              isMatch = isSublist(labelList, matchLabelList);
            } else {
              isMatch = isSublist(matchLabelList, labelList);
            }

            if (isMatch) {
              matches.push(label);
            }
          }

          let bestMatch = "";
          let bestMatchCoverage = 0;
          // If there are matches, find the one with the highest coverage
          if (matches.length > 0) {
            for (let label of matches) {
              let labelList = tokenize(label);

              // Get the coverage of the match label in the group label
              let coverage = longestCommonRun(labelList, matchLabelList);

              if (coverage > bestMatchCoverage) {
                bestMatchCoverage = coverage;
                bestMatch = label;
              }
            }
          }

          // If there is a best match, set the colour to the same as the best match
          if (bestMatch !== "") {
            g.colour = labeledGroupsColours[bestMatch];
          }
        }
      }
    }

    data.groups = [];
    data.groups = updatedGroups;

    plot.call(chart);

    boldGroupSetUp();
  });

  // Populate label type multi select and bind change handler
  let labelTypes = new Set()
  labelTypes.add("all_label_file") // Add labels here so that it's always on top of the dropdown list
  labelTypes.add("psi_blast_match") 
  labelTypes.add("match_short")
  labelTypes.add("match_very_short")
  data.clusters.forEach(cluster => {
    cluster.loci.forEach(locus => {
      locus.genes.forEach(gene => Object.keys(gene.names).forEach(type => labelTypes.add(type)))
    })
  })

  const changeLabel = event => {
    const type = event.target.value
    d3.selectAll("text.geneLabel")
      .each(d => d.label = d.names[type] ? d.names[type] : d.label)
    update({})
  }
  let select = d3.select("#select-label-type")
    .on("change", changeLabel)
  select.selectAll("option")
    .data(labelTypes)
    .join("option")
    .attr("value", d => d)
    .text(d => d.replaceAll("_", " ")) // Modified to change _ to " " to look nicer

  // Figure layout
  d3.select("#input-scale-factor")
    .on("change", function() {update({plot: {scaleFactor: +this.value}})})
  d3.select("#input-cluster-spacing")
    .on("change", function() {update({cluster: {spacing: +this.value}})})
  d3.select("#input-scale-genes")
    .on("change", function() {update({plot: {scaleGenes: d3.select(this).property("checked")}})})
  // Reset Plot Transform -----------------------------------------------------------------------------------------------
  //d3.select("#input-default-display")
  //  .on("click", function() {
  //    document.getElementsByClassName("clusterMapG")[0].setAttribute("transform", "translate(20,50) scale(1.2)");
  //    plot.call(chart);
  //  })
  //---------------------------------------------------------------------------------------------------------------------

  // Cluster
  d3.select("#input-cluster-align-labels")
    .on("change", function() {update({cluster: {alignLabels: d3.select(this).property("checked")}})})
  d3.select("#input-cluster-hide-coords")
    .on("change", function() {update({cluster: {hideLocusCoordinates: d3.select(this).property("checked")}})})
  d3.select("#input-cluster-name-size")
    .on("change", function() {update({cluster: {nameFontSize: +this.value}})})
  d3.select("#input-locus-name-size")
    .on("change", function() {update({cluster: {lociFontSize: +this.value}})})
  d3.select("#input-locus-spacing")
    .on("change", function() {update({locus: {spacing: +this.value}})})

  // Gene polygon shape
  d3.select("#input-body-height")
    .on("change", function() {update({gene: {shape: {bodyHeight: +this.value}}})})
  d3.select("#input-tip-height")
    .on("change", function() {update({gene: {shape: {tipHeight: +this.value}}})})
  d3.select("#input-tip-length")
    .on("change", function() {update({gene: {shape: {tipLength: +this.value}}})})
  d3.select("#input-gene-stroke-colour")
    .on("change", function() {update({gene: {shape: {stroke: this.value}}})})
  d3.select("#input-gene-stroke-width")
    .on("change", function() {update({gene: {shape: {strokeWidth: +this.value}}})})

  // Gene labels
  d3.select("#input-gene-labels")
    .on("change", function() {update({gene: {label: {show: d3.select(this).property("checked")}}})})
  d3.select("#input-label-rotation")
    .on("change", function() {update({gene: {label: {rotation: +this.value}}})})
  d3.select("#input-label-start")
    .on("change", function() {update({gene: {label: {start: +this.value}}})})
  d3.select("#select-label-anchor")
    .on("change", function() {update({gene: {label: {anchor: this.value}}})})
  d3.select("#input-label-size")
    .on("change", function() {update({gene: {label: {fontSize: +this.value}}})})
  d3.select("#input-label-spacing")
    .on("change", function() {update({gene: {label: {spacing: +this.value}}})})
  d3.select("#select-label-position")
    .on("change", function() {update({gene: {label: {position: this.value}}})})
  d3.select("#input-clear-labels")
    .on("click", function() {
      d3.selectAll("text.geneLabel").each(d => d.label = " ");
      document.getElementById("select-label-type").value = "clear";
      plot.call(chart);
    })

  // Scale bar
  d3.select("#input-scalebar-fontsize")
    .on("change", function() {update({scaleBar: {fontSize: +this.value}})})
  d3.select("#input-scalebar-height")
    .on("change", function() {update({scaleBar: {height: +this.value}})})
  d3.select("#input-scalebar-basepair")
    .on("change", function() {update({scaleBar: {basePair: +this.value}})})
  d3.select("#input-scalebar-margin-top")
    .on("change", function() {update({scaleBar: {marginTop: +this.value}})})

  // Colour bar
  d3.select("#input-colourbar-fontsize")
    .on("change", function() {update({colourBar: {fontSize: +this.value}})})
  d3.select("#input-colourbar-height")
    .on("change", function() {update({colourBar: {height: +this.value}})})
  d3.select("#input-colourbar-margin-top")
    .on("change", function() {update({colourBar: {marginTop: +this.value}})})

  // Legend
  d3.select("#input-legend-fontsize")
    .on("change", function() {update({legend: {fontSize: +this.value}})})
  d3.select("#input-legend-entryheight")
    .on("change", function() {update({legend: {entryHeight: +this.value}})})
  d3.select("#input-legend-margin-left")
    .on("change", function() {update({legend: {marginLeft: +this.value}})})

  // Dynamic Gene Labelling -------------------------------------------------------------------------------------------
  // Show/hide dynamic labelling options, if checked, adds "dynamic match" and "dynamic short" to list of label types in dropdown menu
  /**
  d3.select("#input-use-dynamic-labels")
    .on("change", function() {
      let hide = document.getElementById("hide-setting-dynamic")
      if (d3.select(this).property("checked")) {
        hide.setAttribute("style", "")
      } else {
        hide.setAttribute("style", "display: none")
      }
    });
    **/

  // Select file to use as master label file, loop through and assign dynamic match and dynamic short
  d3.select("#input-label-name-path")
    .on("change", function() {

      // Check file size
      let file = document.getElementById("input-label-name-path").files[0];
      if (file && file.size > 500 * 1024 * 1024) { // 500MB
        alert("File is too large. Please upload a file smaller than 500MB.");
        this.value = "";
        file = null;
        return;
      }

      // Display loading div
      document.getElementById('loading').style.display = 'block';

      let existingAccessionsLookup = {};

      for (cluster of data.clusters) {
        for (locus of cluster.loci) {
          for (gene of locus.genes) {
            if ('protein_id' in gene.names) {
              let geneAcc = gene.names['protein_id'];
              if (! (geneAcc in existingAccessionsLookup)) {
                existingAccessionsLookup[geneAcc] = gene;
              } else {
                console.log("Warning: Gene accession '" + geneAcc + "' in the data object multiple times. Retaining only the first occurrence.");
              }
            }
          }
        }
      }

      if (! labelTypes.has("dynamic_match")) {
        // Add dynamic options to dropdown menu
        labelTypes.add("dynamic_match");
        labelTypes.add("dynamic_short");
        labelTypes.add("dynamic_very_short");

        let match = document.createElement("option")
        let short = document.createElement("option")
        let vshort = document.createElement("option")

        match.setAttribute("value", "dynamic_match");
        match.textContent = "dynamic match"

        short.setAttribute("value", "dynamic_short")
        short.textContent = "dynamic short"

        vshort.setAttribute("value", "dynamic_very_short")
        vshort.textContent = "dynamic very short"

        document.getElementById("select-label-type").appendChild(match)
        document.getElementById("select-label-type").appendChild(short)
        document.getElementById("select-label-type").appendChild(vshort)
      }

      // ---------------------------- TODO: SHOULD CHECK IF EACH ACC IN THE LABEL FILE EXISTS IN THE PLOT BEFORE ADDING IT ----------------------------
      // Read the label file into a format with dict lookup (e.g., accession: [label, evalue, [class1, colour1, class2, colour2, ...]])
      // By only keeping accs present in the plot we'll limit the size of the acc dict
      let fr = new FileReader();
      fr.onload=function(){
        let dynamicAccessions = {};
        let baseHeader = "#All_Accessions,Match_Label,Match_E_Value";

        let t = fr.result;
        let lines = t.split(/\r?\n/);

        if (t.slice(0, baseHeader.length) != baseHeader) {
          // File doesn't start with the valid baseHeader, display error and exit
          alert("Error: First line of uploaded file does not begin with \"#All_Accessions,Match_Label,Match_E_Value\"");
          console.log("Error: First line of uploaded file does not begin with \"#All_Accessions,Match_Label,Match_E_Value\"");
          document.getElementById("loading").style.display = 'none';
          return;
        }

        // Header begins with the valid base, add dynamic tier names to list
        let header = lines[0].split(",");
        let tierNames = [];
        for (let i = 3; i < header.length; i++) {
          tierNames.push("Dynamic_" + header[i]);
        }

        for (let i = 1; i < lines.length; i++) { // Start at i = 1 to skip header line
          let splitLine = lines[i].split(',');
          if (splitLine.length < 3) {
            console.log("Warning: Line at index '" + i.toString() + "' has length less than 3 (" + splitLine.length + ") when split on ','. Skipping...");
          } else {
            let accession = splitLine[0];
            let label = splitLine[1];
            let evalue = splitLine[2];

            // Make sure the line length is the same as the header (shouldn't be necessary)
            while (splitLine.length < header.length) {
              console.log("Warning: Line at index '" + i.toString() + "' has less items than the header. Padding with ''...");
              splitLine.push("");
            }

            let lineTiers = splitLine.slice(3,);

            if (accession in dynamicAccessions) {
              if (dynamicAccessions[accession][1] > evalue) { // Duplicate accession (shouldn't really happen), keep one with lower e-value
                // If one doesn't have an e-value, but has tier labels, and the other has an e-value, but no tier labels, move the tier labels over to the new one
                replaceDupe = true;
                let existing = dynamicAccessions[accession];
                let existingTiers = existing[2];

                filledTiers = [];
                for (let i=0; i < lineTiers.length; i++) {
                  if (lineTiers[i]) {
                    // lineTiers[i] is non-empty
                    filledTiers.push(lineTiers[i]);
                  } else if (! lineTiers[i] && existingTiers[i]) {
                    // lineTiers[i] is "" and existingTiers[i] isn't ""
                    filledTiers.push(existingTiers[i]);
                  } else {
                    // Both are empty
                    filledTiers.push("");
                  }
                }

                // Make sure lengths of filledTiers == existingTiers == lineTiers
                if (! (filledTiers.length == existingTiers.length == lineTiers.length)) {
                  console.log("Warning: Length of filledTiers (" + filledTiers.length + "), existingTiers (" + existingTiers.length + "), and lineTiers (" + lineTiers.length + ") do not match. Using lineTiers.");
                  filledTiers = lineTiers;
                }

                dynamicAccessions[accession] = [label, evalue, filledTiers];
              } else {
                // Also check the tiers here
                let existing = dynamicAccessions[accession];
                let existingTiers = existing[2];

                filledTiers = [];
                for (let i=0; i < lineTiers.length; i++) {
                  if (existingTiers[i]) {
                    //existingTiers[i] is non-empty
                    filledTiers.push(existingTiers[i]);
                  } else if (! existingTiers[i] && lineTiers[i]) {
                    // existingTiers[i] is "" and lineTiers[i] isn't ""
                    filledTiers.push(lineTiers[i]);
                  } else {
                    // Both are emtpy
                    filledTiers.push("");
                  }
                }

                // Make sure lengths of filledTiers == existingTiers == lineTiers
                if (! (filledTiers.length == existingTiers.length == lineTiers.length)) {
                  console.log("Warning: Length of filledTiers (" + filledTiers.length + "), existingTiers (" + existingTiers.length + "), and lineTiers (" + lineTiers.length + ") do not match. Using existingTiers.");
                  filledTiers = existingTiers;
                }

                dynamicAccessions[accession] = [existing[0], existing[1], filledTiers];
                continue;
              }
            } else { // New accession, add line
              dynamicAccessions[accession] = [label, evalue, lineTiers];
            }
          }
        }

        // Build a label map for each tier column
        let n_label_dict = [];

        for (let i=0; i < tierNames.length; i++) {
          n_label_dict.push({});
        }

        // Loop through all the genes in the plot, and build the label dicts and colour maps (i.e., parse out the label file data for the genes present in the plot)
        for (accession in existingAccessionsLookup) {
          // Iterate over the keys (accessions) in the plot and check if they're in the label file
          if (accession in dynamicAccessions) {
            // Fill out the label map for each tier column
            for (let i = 0; i < tierNames.length; i++) {
              let tierLabel = dynamicAccessions[accession][2][i];
              let tierColour = "";

              if (tierLabel && tierLabel.trim()) {
                if (tierLabel in n_label_dict[i]) {
                  n_label_dict[i][tierLabel].push(existingAccessionsLookup[accession].uid);
                } else {
                  n_label_dict[i][tierLabel] = [existingAccessionsLookup[accession].uid];
                }
              }
            }
          }
        }

        // Done extracting the specific label file data for the genes in the plot

        // Build the default (PID homology) label and colour maps for each tier
        let n_uid_fn = [];
        let n_default_colour_map = [];
        let n_default_labels = [];
        let n_default_group_order = [];

        for (let i=0; i < tierNames.length; i++) {
          let uid_fn = {};
          let default_colour_map = {};
          let default_labels = {};
          let default_group_order = [];

          for (group of data.groups) {
            // Preserve the order of the groups in the label file (so that re-colouring doesn't shuffle the legend)
            default_group_order.push(group.label);

            // Get UIDs for each gene in the group
            for (geneUID of group.genes) {
              if (group.label in uid_fn) {
                uid_fn[group.label].push(geneUID);
              } else {
                uid_fn[group.label] = [geneUID];
              }
            }

            let strippedColour = null;
            if (group.colour && group.colour.length > 3) { // Check if the colour is RGB
              if (group.colour.includes("rgb(")) {
                strippedColour = group.colour;
                default_colour_map[group.label] = [strippedColour];
              }
            }

            if (!strippedColour || !strippedColour.trim()) { // If colour isn't RGB, check if it's a hex colour
              strippedColour = "#" + group.colour.replace(/[^0-9a-fA-F]/g, "");
              if (strippedColour != "#" && strippedColour != null && strippedColour != " ") {
                default_colour_map[group.label] = [strippedColour];
              } else {
                default_colour_map[group.label] = [""];
              }
            }

            default_labels[group.label] = group.label;
          }
          
          // Push the column data to the n storage arrays
          n_uid_fn.push(uid_fn);
          n_default_colour_map.push(default_colour_map);
          n_default_labels.push(default_labels);
          n_default_group_order.push(default_group_order);
        }

        // Now we have all the data we need to build the new groups (homology extended and non-homologous subclass)
        // 1 - If a gene from the label file is in a homology group, rename the homology group to the label (e.g., Group 3 --> "Abi2 (dynamic homologous)")
        // 2 - If another gene with a different label from the label file is in the same homology group, prepend it to the label (e.g., Abi2 (dynamic homologous) --> Abi3 | Abi2 (dynamic homologous))
        // 3 - If a gene from the label file isn't in a homology group, create a new group for it in the fn and colour it black
        //   - If the label is already in the fn, add the gene to it (e.g., if 3 already happened)

        let keys = Object.keys(n_uid_fn[0]);
        for (let i = 0; i < n_label_dict.length; i++) {
          let labelDict = n_label_dict[i];
          let colourMap = n_default_colour_map[i];
          let names = n_default_labels[i];
          let defaultFN = n_uid_fn[i];

          let uid_to_labels = {};
          let uid_to_group = {};

          // Flip the keys and values in the two dicts
          for (let newLabel in labelDict) {
            let groupUIDs = labelDict[newLabel];
            for (let uid of groupUIDs) {
              uid_to_labels[uid] = newLabel;
            }
          }

          for (let groupLabel in defaultFN) {
            let groupUIDs = defaultFN[groupLabel];
            for (let uid of groupUIDs) {
              uid_to_group[uid] = groupLabel;
            }
          }

          // Get a list of all the UIDs in the label dict
          let labelUIDs = [];
          for (let label in labelDict) {
            for (let uid of labelDict[label]) {
              labelUIDs.push(uid);
            }
          }

          let groupUIDs = [];
          for (let groupLabel in defaultFN) {
            for (let uid of defaultFN[groupLabel]) {
              groupUIDs.push(uid);
            }
          }

          // Get the intersection of the two lists
          let intersection = labelUIDs.filter(uid => groupUIDs.includes(uid));
          let labelRest = labelUIDs.filter(uid => !groupUIDs.includes(uid));

          for (let uid of intersection) {
            let label = uid_to_labels[uid];
            let group = uid_to_group[uid];

            if (names[group] == group) {
              // Group name hasn't been updated, overwrite it
              names[group] = label + "(Dynamic Homologous " + group.split(": ")[0] + ")";
            } else {
              // Group name has been updated, check if it has the label already
              if (!(names[group].includes(label + " (Dynamic Homologous ") || names[group].includes(label + " | "))) {
                // Group name has been updated with a different label, add this one to it
                names[group] = label + " | " + names[group];
              }
            }
          }

          let remove = []; // List of entries to remove from fn (those that were renamed)

          // Update colour map with new labels
          for (let group of keys) {
            if (names[group] == group) {
              // Group name is not updated, colour it grey
              colourMap[group] = ["#808080"];
            } else {
              defaultFN[names[group]] = defaultFN[group];
              remove.push(group);
              // Update colour map with new label
              colourMap[names[group]] = colourMap[group];
            }
          }

          // Remove the old labels from the maps
          for (let group of remove) {
            delete defaultFN[group];
            delete colourMap[group];
          }

          // Add the rest of the label genes to the fn
          for (uid of labelRest) {
            label = uid_to_labels[uid] + " (Dynamic)";

            if (label in defaultFN) {
              // New group already in fn, add the uid to it
              defaultFN[label].push(uid);
            } else {
              // Add label to fn, set colour to black
              defaultFN[label] = [uid];
              colourMap[label] = ["#000000"];

              // Add the label to the order list
              n_default_group_order[i].push(label);
            }
          }

          // Homology (group) genes should already be in there
        }

        // Go through the FNs/Colour Maps and build the groups
        for (let i=0; i < n_uid_fn.length; i++) {
          fn = n_uid_fn[i];
          colourMap = n_default_colour_map[i];
          groups = [];
          order = n_default_group_order[i];
          names = n_default_labels[i];

          for (let groupLabel of order) {
            let labeledGroup = false; // Flag to disable GroupNum toggling
            if (groupLabel in names) {
              if (names[groupLabel] != groupLabel || names[groupLabel].includes("Homologous")) {
                labeledGroup = true; // Group name updated
              }
              groupLabel = names[groupLabel]
            } else {
              labeledGroup = true; // New labeled group
            }
            
            if (groupLabel in colourMap && colourMap[groupLabel].length > 1) {
              // If the group label has more than one colour, set it to the first one
              alert("Warning: Group label " + groupLabel + " has multiple colours assigned to it. Using the first one in the label file.");
            } else if (!(groupLabel in colourMap)) {
              alert(groupLabel + " not in colour map");
              colourMap[groupLabel] = ["#808080"];
            }

            colourMap[groupLabel] = colourMap[groupLabel][0];

            newGroup = {
              "uid": "dynamic-group" + i + "-" + groupLabel.split(": ")[0].trim().replace(" ", ""),
              "label": groupLabel,
              "genes": fn[groupLabel],
              "hidden": true,
              "colour": colourMap[groupLabel],
              "labeled": labeledGroup,
              "baseLabel": groupLabel.split(": ")[0].trim()
            };

            // Add the group to the groups array
            groups.push(newGroup);
          }

          data.classes.push({"name": tierNames[i], "groups": groups}); // Add the new classes to the data object
        }



        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


        // Assign the labels to the genes
        let matchDict = dynamicAccessions;
        // let newFlag = isNew

        let classDicts = [];

        for (let i=0; i < n_label_dict.length; i++) {
          // Flip the keys and values in the label dict
          let labelDict = n_label_dict[i];
          let classDict = {};
          for (let label in labelDict) {
            let labelUIDs = labelDict[label];
            for (let uid of labelUIDs) {
              classDict[uid] = label;
            }
          }
          classDicts.push(classDict);
        }

        // Build the new dynamic label names
        let labelNames = [];
        for (let i=0; i<classDicts.length; i++) {
          let labelName = tierNames[i];
          labelName = labelName.replace("Dynamic_", "");
          labelName += "_(dynamic_custom_" + i + ")";
          labelNames.push(labelName);
        }

        // Add the new label names to the dropdown
        for (let name of labelNames) {
          let newLabel = document.createElement("option");
          newLabel.setAttribute("value", name);
          newLabel.textContent = name.replaceAll("_", " ");
          document.getElementById("select-label-type").appendChild(newLabel);
        }

        for (cluster of data.clusters) {
          for (locus of cluster.loci) {
            for (gene of locus.genes) {
              let flag = false;
              let geneUID = gene.uid;

              // Label file is in the new format, add text class labels
              for (let i=0; i<classDicts.length; i++) {
                if (geneUID in classDicts[i]) {
                  gene.names[labelNames[i]] = classDicts[i][geneUID];
                } else {
                  gene.names[labelNames[i]] = "";
                }
              }

              // Add the normal dynamic text labels
              if ("protein_id" in gene.names && gene.names["protein_id"] in matchDict) {
                const pId = gene.names["protein_id"]
                gene.names["dynamic_match"] = pId + " (" + matchDict[pId][0] + ", " + matchDict[pId][1] + ")";
                gene.names["dynamic_short"] = matchDict[pId][0] + " (" + matchDict[pId][1] + ")";
                gene.names["dynamic_very_short"] = matchDict[pId][0]
              } else { flag = true }

              if (flag) {
                gene.names["dynamic_match"] = "";
                gene.names["dynamic_short"] = "";
                gene.names["dynamic_very_short"] = "";
              }
            } 
          }
        }

        // Add the new class labels to the dropdown
        for (let i = 0; i < tierNames.length; i++) {
          let newClassName = tierNames[i];
          let newOption = document.createElement("option");
          newOption.setAttribute("value", newClassName);
          newOption.textContent = newClassName;
          document.getElementById("select-label-class").appendChild(newOption);

        }

        d3.selectAll("text.geneLabel").each(d => d.label = " ");
        let e = document.getElementById("select-label-type");
        let type = e.options[e.selectedIndex].value;
        d3.selectAll("text.geneLabel").each(d => d.label = d.names[type] ? d.names[type] : d.label)
        plot.call(chart);

        alert("Dynamic labels loaded.");
        document.getElementById("loading").style.display = "none"; // Hide loading div
      }


    fr.readAsText(document.getElementById("input-label-name-path").files[0]);

  });
  //-------------------------------------------------------------------------------------------------------------------

  // MODIFY GROUP LABELS ------------------------------------------------------------------------------------------------
  d3.select("#input-group-label-show")
  .on("change", function() {
    // Need to check if group.groupNum exists!
    if (d3.select(this).property("checked")) {
      for (let group of data.groups) {
        if (!group.labeled) {
          if (group.label.includes("Group "+group.groupNum)) {}
          else {(group.label != "" && group.label != " ")
          ? group.label = "Group "+group.groupNum+": "+group.label
          : group.label = "Group "+group.groupNum}
        }
      }
    } else {
      for (let group of data.groups) {
        if (!group.labeled) {
          if (group.label.includes("Group "+group.groupNum)) {
            if (group.label.indexOf(": ") != -1) {
              group.label = group.label.split(": ").slice(1,).join(": ")
            } else { group.label = "" }
          }
        }
      }
    }

    plot.call(chart);
  });
  d3.select("#input-group-label-length")
    .on("change", function() {
      updateGroupProductsDisplay(this.value);

      plot.call(chart);
    });
  //---------------------------------------------------------------------------------------------------------------------


  // Links
  d3.select("#input-link-show")
    .on("change", function() {update({link: {show: d3.select(this).property("checked")}})})
  d3.select("#input-link-best-only")
    .on("change", function() {update({link: {bestOnly: d3.select(this).property("checked")}})})
  d3.select("#input-link-as-line")
    .on("change", function() {update({link: {asLine: d3.select(this).property("checked")}})})
  d3.select("#input-link-straight")
    .on("change", function() {update({link: {straight: d3.select(this).property("checked")}})})
  d3.select("#input-link-group-colour")
    .on("change", function() {update({link: {groupColour: d3.select(this).property("checked")}})})
  d3.select("#input-link-stroke-width")
    .on("change", function() {update({link: {strokeWidth: +this.value}})})
  d3.select("#input-link-threshold")
    .on("change", function() {update({link: {threshold: +this.value}})})
  d3.select("#input-link-fontsize")
    .on("change", function() {update({link: {label: {fontSize: +this.value}}})})
  d3.select("#input-link-label-pos")
    .on("change", function() {update({link: {label: {position: +this.value}}})})
  d3.select("#input-link-label-show")
    .on("change", function() {update({link: {label: {show: d3.select(this).property("checked")}}})})
  d3.select("#input-link-labelbg-show")
    .on("change", function() {update({link: {label: {background: d3.select(this).property("checked")}}})})



  //---------------------------------------------------------------------------------------------------------------------
  // Block of custom code to change the default plot values, flip backwards loci, and set the anchor gene

  
  // Changing the clustermap default config for showing gene labels doesn't work, so have to set it here
  d3.selectAll("text.geneLabel").each(d => d.label = d.names["all_label_file"] ? d.names["all_label_file"] : d.names["clear"]); // Default label all label file if present, else clear (" ")
  d3.select("#plot").call(ClusterMap.ClusterMap()
    .config({
      gene: {
        label: {
          show: true,
        }
      },
    }))

    // Simulate mouse events to flip loci and set the anchor gene
  const dblClickEvent = new MouseEvent('dblclick', {
    'view': window,
    'bubbles': true,
    'cancelable': true
  });

  const clickEvent = new MouseEvent('click', {
    'view': window,
    'bubbles': true,
    'cancelable': true
  });

  for (let locus_id of flipLoci) {
    document.getElementById(locus_id).dispatchEvent(dblClickEvent)
  }

  for (anchorGeneUid of anchorGenes) {
    document.getElementById(anchorGeneUid).getElementsByTagName("polygon")[0].dispatchEvent(clickEvent);
  }

  setGroupAttributes();
  boldGroupSetUp();
  updateGroupProductsDisplay(2);

  // End of custom block
  //---------------------------------------------------------------------------------------------------------------------
}


if (typeof data === 'undefined') {
  const data = d3.json("data.json").then(plot)
} else {
  plot(data)
}

//-----------------------------------------------------------------------------------------------------------------------
// More Custom

// Modify title element (name on tab) to display which file is open
let filename = window.location.pathname.split("/").pop().split(".")[0];
document.getElementsByTagName("title")[0].textContent = filename.replace(/_/g, " ") + " (clinker)";


// Add event listeners to modify cluster labels on mouse hover (enter/exit) to display the full name or shortened name
const targets = document.getElementsByClassName("clusterText");

for (target of targets) {
    target.setAttribute("data-original-text", target.textContent);

    const splitName = target.getAttribute("data-original-text").split(" ");
    let shortName = "";
    let shortLen;
    if (splitName[0] == "TPA") {
        if (4 > splitName.length) {
            shortLen = splitName.length;
        } else {
            shortLen = 4;
        }
    } else { shortLen = (splitName.length < 2) ? splitName.length : 2; }

    for ( let i = 0; i < shortLen; i++ ) {
        shortName += splitName[i];
        shortName += " ";
    }
    shortName = shortName.slice(0,-1)+"...";

    target.setAttribute("data-short-text", shortName);

    target.textContent = target.getAttribute("data-short-text");

    target.addEventListener("mouseenter", function() {
        this.textContent = this.getAttribute("data-original-text");
    });

    target.addEventListener("mouseleave", function() {
        this.textContent = this.getAttribute("data-short-text");
    });
}

function boldGroup(g) {
  // Function to bold-ify the borders of genes within a group
  g.style.fontWeight = 'bold';
  let gUid = g.getElementsByTagName("circle")[0].getAttribute("class");
  // For some reason, circle class isn't the same as the polygon class
  //   e.g., "group-dynamic-group0-Group-0" vs "group-dynamic-group0-Group0"
  for (e of document.getElementsByClassName(gUid)) {
    // alert(e.tagName);
    if (e.tagName == "polygon") {
      e.setAttribute("style", "stroke: black; stroke-width: 5")
    }
  }
}

function unboldGroup(g) {
  // Function to unbold-ify the borders of genes within a group
  let gUid = g.getElementsByTagName("circle")[0].getAttribute("class");
  g.style.fontWeight = 'normal';
  for (e of document.getElementsByClassName(gUid)) {
    if (e.tagName == "polygon") {
      e.setAttribute("style", "stroke: black; stroke-width: 1")
    }
  }
}

function checkBold(g) {
  // Function to check if a group is bold-ified
  // If any of the genes in the group are bolded, return true
  // If the group name itself is bolded, return true
  let gUid = g.getElementsByTagName("circle")[0].getAttribute("class");
  for (e of document.getElementsByClassName(gUid)) {
    if (e.tagName == "polygon") {
      if (e.getAttribute("style").includes("stroke-width: 5")) {
        return true;
      }
    }
  }
  return g.style.fontWeight == 'bold';
}

function boldGroupSetUp() {
  // Add event listeners to thicken borders of genes that are members of hovered group
  let groupTargets = document.getElementsByClassName("legend")[0].getElementsByClassName("element");

  for (g of groupTargets) {
    g.boldToggled = false; // Set a flag to check if the group was bolded via dblclick

    g.addEventListener("mouseenter", function() {
      if (!checkBold(this) && !this.boldToggled) {
        boldGroup(this);
      }
    })

    g.addEventListener("mouseleave", function() {
      if (checkBold(this) && !this.boldToggled) {
        unboldGroup(this);
      }
    })

    g.addEventListener("click", function(e) {
      if (e.shiftKey) {
        if (this.boldToggled) {
          this.boldToggled = false;
          unboldGroup(this);
        } else {
          this.boldToggled = true;
          boldGroup(this);
        }
      }
    })
  }
}

function resetBolding() {
  // Reset bolding of all groups (called when group type is changed)
  const groupTargets = document.getElementsByClassName("legend")[0].getElementsByClassName("element");

  for (g of groupTargets) {
    if (checkBold(g)) {
      unboldGroup(g);
    }
  }
}

// Add group event listeners to groups on page load
// boldGroupSetUp();
// document.addEventListener("DOMContentLoaded", boldGroupSetUp);

function setGroupNums(groups=data.groups) {
  // Function to set the group.groupNum attribute for each group (just use the index of the group in the array)
  let groupNum = 0;
  for (let group of groups) {
    // Check if group.groupNum exists
    if (!group.labeled) {
      group.groupNum = groupNum;
    }
    groupNum++;
  }
}

function setGroupProducts(groups=data.groups) {
  // Function to populate the group.products attribute for each group
  for (let group of groups) {
    let productFreqs = {};
    for (let geneUID of group.genes) {
      // Need to get the actual gene object from the data.clusters array
      for (let cluster of data.clusters) {
        for (let locus of cluster.loci) {
          for (let gene of locus.genes) {
            if (gene.uid == geneUID && "product" in gene.names) {
              let product = gene.names["product"];

              let disallowedProducts = ["hypothetical protein", "uncharacterized protein"]

              if (!(product.toLowerCase() in disallowedProducts) && product.trim() != "") {
                if (product in productFreqs) {
                  productFreqs[product]++;
                } else {
                  productFreqs[product] = 1;
                }
              }

              // Check if the gene is a PSI-BLAST match
              if ("psi_blast_match" in gene.names && gene.names["psi_blast_match"] != null && gene.names["psi_blast_match"].trim() != "") {
                if (product in productFreqs) {
                  productFreqs[gene.names["match_very_short"]] = 999;
                } else {
                  productFreqs[gene.names["match_very_short"]] = 999;
                }
              }
            }
          }
        }
      }
    }

    let sortedProducts = Object.entries(productFreqs).sort((a, b) => b[1] - a[1]);
    let productList = [];
    for (let i = 0; i < sortedProducts.length; i++) {
      productList.push(sortedProducts[i][0]);
    }

    // If the group contains a PSI-BLAST match, prepend the Query name to the product list


    group.products = productList;
  }
}

function updateGroupProductsDisplay(n) {
  newLen = n;

  for (let group of data.groups) {
    // Check if group.products exists
    if (!group.products) {
      group.products = [];
    }

    let newArr = []
    for (let i = 0; (i < newLen && i < group.products.length); i++) {
      newArr.push(group.products[i])
    }
    //alert([newArr, newLen, group.products.length])
    if (d3.select("#input-group-label-show").property("checked")) {
      group.label = group.baseLabel + ": "+newArr.join(" | ");
    } else {
      if (!group.labeled) {
        group.label = newArr.join(" | ");
      } else {
        // If it's a label group, add the products to the baseLabel
        group.label = group.baseLabel + ": " + newArr.join(" | ");
      }
    }
    if (group.label[group.label.length -2] == ":") { group.label = group.label.slice(0, group.label.length-2)}
  }
}

function setGroupAttributes() {
  // Function to add the new group attributes to each group
  if (data.classes !== undefined) {
    for (let tier of data.classes) {
      setGroupProducts(tier.groups);
      for (let group of tier.groups) {
        let baseLabel = group.label.split(": ")[0].trim();
        group.baseLabel = baseLabel;

        if (baseLabel.startsWith("Group ")) {
          group.labeled = false;
        } else {
          group.labeled = true;
        }
      }
      setGroupNums(tier.groups);
    }
  }


  setGroupProducts(data.groups);
  for (let group of data.groups) {
    let baseLabel = group.label.split(": ")[0].trim();
    group.baseLabel = baseLabel;

    if (baseLabel.startsWith("Group ")) {
      group.labeled = false;
    } else {
      group.labeled = true;
    }
  }
  setGroupNums(data.groups);
}


// End
//-----------------------------------------------------------------------------------------------------------------------


