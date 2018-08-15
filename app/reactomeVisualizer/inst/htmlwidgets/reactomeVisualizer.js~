HTMLWidgets.widget({

  name: 'reactomeVisualizer',

  type: 'output',

    factory: function(el, width, height) {

    // TODO: define shared variables for this instance

	return {
	//Creating the Reactome pathways overview widget
	    onReactomeFireworksReady: function(){
	    //This function is automatically called when the widget code is ready to be used
		var fireworks = Reactome.Fireworks.create({
		    "placeHolder" : el,
		    "width" : width,
		    "height" : height
		});
		
		//Adding different listeners
		
		fireworks.onFireworksLoaded(function (loaded) {
		    console.info("Loaded ", loaded);
		});
		
		fireworks.onNodeHovered(function (hovered){
		    console.info("Hovered ", hovered);
		});
	    
		fireworks.onNodeSelected(function (selected){
		    console.info("Selected ", selected);
		});
		fireworks.highlightNode(x.pathId)
	    };
	}
    };
  }
});
