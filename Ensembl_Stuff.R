ensg_map_obj <- readRDS("Ensembl_name_maps.rds")
ensg_name_map <- ensg_map_obj[["ensg2symbol"]]
ensg2musg <- ensg_map_obj[["ensg2musg"]]
ensg2rat <- ensg_map_obj[["ensg2rat"]]
musg2rat <- ensg_map_obj[["musg2rat"]]
musg_name_map <- ensg_map_obj[["musg2symbol"]]
rat_name_map <- ensg_map_obj[["rat2symbol"]]

h_dups = unique(ensg2musg[(duplicated(ensg2musg[,1])),1])
m_dups = unique(ensg2musg[(duplicated(ensg2musg[,2])),2])
ensg2musg_one2one <- ensg2musg[!(ensg2musg[,1] %in% h_dups) & !(ensg2musg[,2] %in% m_dups),]

map_symbol_ensg <- function(genes, is.org=c("Hsap","Mmus","Rat"), is.name=c("symbol","ensg")) {
	if (is.org[1] == "Hsap") {
		if(is.name[1]=="symbol") {
			new = as.character(ensg_name_map[match(genes, ensg_name_map[,2]),1])
		} else if (is.name[1] =="ensg") {
			new = as.character(ensg_name_map[match(genes, ensg_name_map[,1]),2])
		} else {
			stop("Unrecognized name type")
		}
	} else if (is.org[1] == "Mmus") {
		if(is.name[1]=="symbol") {
			new = as.character(musg_name_map[match(genes, musg_name_map[,2]),1])
		} else if (is.name[1] =="ensg") {
			new = as.character(musg_name_map[match(genes, musg_name_map[,1]),2])
		} else {
			stop("Unrecognized name type")
		}
	} else if (is.org[1] == "Rat") {
		if(is.name[1]=="symbol") {
			new = as.character(rat_name_map[match(genes, rat_name_map[,2]),1])
		} else if (is.name[1] =="ensg") {
			new = as.character(rat_name_map[match(genes, rat_name_map[,1]),2])
		} else {
			stop("Unrecognized name type")
		}
	} else {
		stop("Unrecognized organism");
	}
#	new[is.na(new)] = as.character(genes[is.na(new)])
	new[is.na(new)] = ""
	return(new);
}

map_Hsap_Mmus <- function(genes, is.org=c("Hsap","Mmus"), one2one=FALSE) {
	if (one2one) {
		local_map <- ensg2musg_one2one
	} else {
		local_map <- ensg2musg
	}
	if (is.org[1] == "Hsap") {
		new = as.character(local_map[match(genes, local_map[,1]),2])
	} else if (is.org[1] == "Mmus") {
		new = as.character(local_map[match(genes, local_map[,2]),1])
	} else {
		stop("Unrecognized organism");
	}
#	new[is.na(new)] = as.character(genes[is.na(new)])
	new[is.na(new)] = ""
	return(new);
}

map_ensg_across_org <- function(genes, is.org=c("Hsap", "Mmus", "Rat"), out.org=c("Hsap", "Mmus", "Rat"), one2one=FALSE) {
	if (is.org == "Hsap") {
		if (out.org == "Mmus") {
			cross_spp_map <- ensg2musg
		} else if (out.org == "Rat") {
			cross_spp_map <- ensg2rat
		} else {
			stop("Unrecognized organism")
		}
	} else if (is.org == "Mmus") {
		if (out.org == "Hsap") {
			cross_spp_map <- cbind(as.character(ensg2musg[,2]), as.character(ensg2musg[,1]))
		} else if (out.org == "Rat") {
			cross_spp_map <- musg2rat;
		} else {
                        stop("Unrecognized organism")
                }
	} else if (is.org == "Rat") {
		 if (out.org == "Mmus") {
                        cross_spp_map <- cbind(as.character(musg2rat[,2]), as.character(musg2rat[,1]))
                } else if (out.org == "Hsap") {
                        cross_spp_map <- cbind(as.character(ensg2rat[,2]), as.character(ensg2rat[,1]))
                } else {
                        stop("Unrecognized organism")
                }
	} else {
		stop("Unrecognized organism")
	}
	if (one2one) {
		first_dup <- unique(cross_spp_map[duplicated(cross_spp_map[,1]),1])
		second_dup <- unique(cross_spp_map[duplicated(cross_spp_map[,2]),2])
		cross_spp_map <- cross_spp_map[!(cross_spp_map[,1] %in% first_dup) &
						!(cross_spp_map[,2] %in% second_dup), ]
	}

	new = as.character(cross_spp_map[match(genes, cross_spp_map[,1]),2])
	new[is.na(new)] = ""
	return(new);
}

General_Map <- function(genes, in.org=c("Hsap","Mmus","Rat"), in.name=c("symbol","ensg"), out.org=c("Hsap","Mmus", "Rat"), out.name=c("symbol","ensg")) {
	if (in.org == out.org & in.name == out.name) {
		# No change
		return(genes)
	}
	if (in.org == out.org) {
		# change names not spp
		return(map_symbol_ensg(genes, is.org=in.org, is.name=in.name))

	} else {
		if (in.name == "symbol") {
			tmp <- map_symbol_ensg(genes, is.org=in.org, is.name=in.name)
		} else {
			tmp <- genes
		}
		tmp <- map_ensg_across_org(tmp, is.org=in.org, out.org=out.org)
		if (out.name =="symbol") {
			out <- map_symbol_ensg(tmp, is.org=out.org, is.name="ensg")
		} else {
			out <- tmp
		}
		return(out)	
	}
}

M_symbol_2_H_symbol <- function(genes) {
        tmp <- map_symbol_ensg(genes, is.org="Mmus", is.name="symbol")
        tmp <- map_Hsap_Mmus(tmp, is.org="Mmus")
        tmp <- map_symbol_ensg(tmp, is.org="Hsap", is.name="ensg")
        return(tmp)
}

H_symbol_2_M_symbol <- function(genes) {
        tmp <- map_symbol_ensg(genes, is.org="Hsap", is.name="symbol")
        tmp <- map_Hsap_Mmus(tmp, is.org="Hsap")
        tmp <- map_symbol_ensg(tmp, is.org="Mmus", is.name="ensg")
        return(tmp)
}

