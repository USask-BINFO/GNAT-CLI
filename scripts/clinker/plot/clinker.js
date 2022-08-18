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


for (g in data.groups) {
  const group = data.groups[g];
  let psiFlag = false;
  let groupProducts = {};
  let newLabel = group.label + ": ";

  for (geneUid of group.genes) {
    if (Object.keys(matchGenes).includes(geneUid)) {
      if (!newLabel.includes(matchGenes[geneUid])) {
        newLabel += matchGenes[geneUid];
        psiFlag = true;
      }
    }
    if (geneUidMap[geneUid][0] != "undefined") {
      groupProducts[geneUidMap[geneUid][0]] == null ? groupProducts[geneUidMap[geneUid][0]] = 1 : groupProducts[geneUidMap[geneUid][0]]++;
    }
    else { alert([geneUidMap[geneUid], geneUid]) }
  }


  let keys = Object.keys(groupProducts);
  let noUndefined = [];
  for (key of keys) {
    if (key != "undefined") { noUndefined.push(key); }
    else { noUndefined.push("pseudogene"); }
  }
  keys = noUndefined;
  let newProducts = {};
  for (key of keys) {
    if (key == "pseudogene") {
      newProducts[key] = groupProducts["undefined"]
    }
    newProducts[key] = groupProducts[key];
  }
  groupProducts = newProducts;

  let first, second;

  if (psiFlag) {
    first = newLabel.slice(newLabel.indexOf(":")+2);
    second = "";
    keys = [];
    for (acr of first.split(",")) {
      acr != "" ? keys.push(acr) : {};
      }
    // data.groups[g].label = newLabel.slice(0, newLabel.length - 1);
  }
  else {
    if (keys.length == 1) {  }
    else {
      const lenGroup = group.genes.length;
      const numKeys = keys.length;

      let sortable = [];
      for (let key in groupProducts) {
        sortable.push([key, groupProducts[key]]);
      }

      sortable.sort(function(a,b) {
        return a[1] - b[1];
      });

      keys = sortable.reverse();
      for (i in keys) {
        keys[i] = keys[i][0];
      }
    }

    if (keys.length > 1) {
        first = keys[0];
        second = keys[1];
    } else if (keys.length > 0) {
        first = keys[0];
        second = "";
    } else {
        first = "";
        second = "";
    }

  }

  if (second != "") {
    data.groups[g].label = group.label + ":" + " " + first + " | " + second;
  } else {
    if (first[first.length - 1] == ",") { first = first.slice(0, first.length-1) }
    data.groups[g].label = group.label + ":" + " " + first;
  }
  data.groups[g].products = keys;

  let gN = group.label.split(" ")[1];
  if (gN[gN.length-1] == ":") { gN = gN.slice(0, gN.length-1) }

  data.groups[g].groupNum = gN;

}

//-----------------------------------------------------------------------------------------------------------------------


function plot(data) {
  const div = d3.select("#plot")
  const chart = ClusterMap.ClusterMap()
    .config({
      scaleFactor: 80,  // Modified, used to be 30
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

  // Populate label type multi select and bind change handler
  let labelTypes = new Set()
  labelTypes.add("psi_blast_match") // Add psi_blast_match here so that it's always on top of the dropdown list
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
      if (! labelTypes.has("dynamic_match")) {

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

      let fr = new FileReader();
      fr.onload=function(){
        let matchDict = {};
        let t = fr.result;
        let lines = t.split(/\r?\n/);
        for (line of lines) {
          const splitLine = line.split("\t");
          if (splitLine.length == 3) {
            matchDict[splitLine[0]] = [splitLine[1], splitLine[2]]
          }
        }
        for (cluster of data.clusters) {
          for (locus of cluster.loci) {
            for (gene of locus.genes) {
              let flag = false;
              if ("protein_id" in gene.names) {
                if (gene.names["protein_id"] in matchDict) {
                  const pId = gene.names["protein_id"]
                  gene.names["dynamic_match"] = pId + " (" + matchDict[pId][0] + ", " + matchDict[pId][1] + ")";
                  gene.names["dynamic_short"] = matchDict[pId][0] + " (" + matchDict[pId][1] + ")";
                  gene.names["dynamic_very_short"] = matchDict[pId][0]
                } else { flag = true }
              } else { flag = true }
              if (flag) {
                gene.names["dynamic_match"] = "";
                gene.names["dynamic_short"] = "";
                gene.names["dynamic_very_short"] = "";
              }
            }
          }
        }
        d3.selectAll("text.geneLabel").each(d => d.label = " ");
        let e = document.getElementById("select-label-type");
        let type = e.options[e.selectedIndex].value;
        d3.selectAll("text.geneLabel").each(d => d.label = d.names[type] ? d.names[type] : d.label)
        plot.call(chart);
      }
      fr.readAsText(document.getElementById("input-label-name-path").files[0]);
    })
  //-------------------------------------------------------------------------------------------------------------------

  // MODIFY GROUP LABELS ------------------------------------------------------------------------------------------------
  d3.select("#input-group-label-show")
  .on("change", function() {
    if (d3.select(this).property("checked")) {
      for (let group of data.groups) {
        if (group.label.includes("Group "+group.groupNum)) {}
        else {(group.label != "" && group.label != " ")
        ? group.label = "Group "+group.groupNum+": "+group.label
        : group.label = "Group "+group.groupNum}
      }
    } else {
      for (let group of data.groups) {
        if (group.label.includes("Group "+group.groupNum)) {
          if (group.label.indexOf(": ") != -1) {
            group.label = group.label.split(": ").slice(1,).join(": ")
          } else { group.label = "" }
        }
      }
    }

    plot.call(chart);
  });
  d3.select("#input-group-label-length")
    .on("change", function() {
      newLen = this.value
      for (let group of data.groups) {
        let newArr = []
        for (let i = 0; (i < newLen && i < group.products.length); i++) {
          newArr.push(group.products[i])
        }
        //alert([newArr, newLen, group.products.length])
        if (d3.select("#input-group-label-show").property("checked")) {
          group.label = "Group "+group.groupNum+": "+newArr.join(" | ");
        } else {
          group.label = newArr.join(" | ");
        }
        if (group.label[group.label.length -2] == ":") { group.label = group.label.slice(0, group.label.length-2)}
      }

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
  function set_defaults() {
  update({gene: {label: {show: true}}}); // Gene labels on by default
  d3.selectAll("text.geneLabel").each(d => d.label = d.names["psi_blast_match"] ? d.names["psi_blast_match"] : d.names["clear"]); update({}); // Default label psi_blast_match if present, else clear (" ")
  update({link: {show: false}}) // Default show links off
  update({plot: {scaleFactor: 80}}) // Default scale 80
  update({gene: {label: {rotation: 45}}}) // Default gene label angle = 45 degrees
  }
  set_defaults()

    // Simulate mouse events to flip loci and set the anchor gene
  const event = new MouseEvent('dblclick', {
    'view': window,
    'bubbles': true,
    'cancelable': true
  });

  const event2 = new MouseEvent('click', {
    'view': window,
    'bubbles': true,
    'cancelable': true
  });

  for (let locus_id of flipLoci) {
    document.getElementById(locus_id).dispatchEvent(event)
  }

  for (anchorGeneUid of anchorGenes) {
    document.getElementById(anchorGeneUid).getElementsByTagName("polygon")[0].dispatchEvent(event2);
  }
  //document.getElementById(anchorGeneId).getElementsByTagName("polygon")[0].dispatchEvent(event2);

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

function updatePlot() {
  // Not used, just in case I need to change something in data and update the plot to reflect the change
  d3.select("#plot").call(ClusterMap.ClusterMap()
  .config({
    scaleFactor: 80,  // Modified, used to be 30
    cluster: {
      spacing: 50,
      alignLabels: true,
    },
    gene: {
      label: {
        show: false,
      }
    },
  }))
}

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

// Add event listeners to thicken borders of genes that are members of hovered group
const groupTargets = document.getElementsByClassName("legend")[0].getElementsByClassName("element");

for (g of groupTargets) {
  g.addEventListener("mouseenter", function() {
    this.style.fontWeight = 'bold';
    let gUid = this.getElementsByTagName("circle")[0].getAttribute("class");
    for (e of document.getElementsByClassName(gUid)) {
      if (e.tagName == "polygon") {
        e.setAttribute("style", "stroke: black; stroke-width: 5")
      }
    }
  })

  g.addEventListener("mouseleave", function() {
    let gUid = this.getElementsByTagName("circle")[0].getAttribute("class");
    this.style.fontWeight = 'normal';
    for (e of document.getElementsByClassName(gUid)) {
      if (e.tagName == "polygon") {
        e.setAttribute("style", "stroke: black; stroke-width: 1")
      }
    }
  })
}



// End
//-----------------------------------------------------------------------------------------------------------------------


