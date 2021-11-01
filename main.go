package main

import (
	"fmt"
	"gifhelper"
	"os"
)

func main() {

	command := os.Args[1]

	g := make([]Galaxy, 0)
	var width float64

	scalingFactor := 1e11 // a scaling factor is needed to inflate size of stars when drawn because galaxies are very sparse

	switch command {
	case "jupiter":
		width = 1e22
		g0 := CreateJupiterSystem(width)
		g = append(g, g0)

		scalingFactor = 1.5e12
	case "galaxy":
		width = 1.0e22
		g0 := InitializeGalaxy(500, 4e21, width/2, width/2)
		g = append(g, g0)

		scalingFactor = 1e10
	case "collide":
		scalingFactor = 1.8e10
		width = 1.0e23
		g0 := InitializeGalaxy(500, 4e21, 7e22, 2e22)
		g1 := InitializeGalaxy(500, 4e21, 3e22, 7e22)
		g0.Push(-1e15, 1e15)
		g1.Push(1e15, -1e15)
		g = append(g, g0, g1)
	}

	// you probably want to apply a "push" function at this point to these galaxies to move
	// them toward each other to collide.
	// be careful: if you push them too fast, they'll just fly through each other.
	// too slow and the black holes at the center collide and hilarity ensues.

	initialUniverse := InitializeUniverse(g, width)

	// now evolve the universe: feel free to adjust the following parameters.
	// numGens := 500000
	numGens := 5000
	time := 2e14
	theta := 0.5

	timePoints := BarnesHut(initialUniverse, numGens, time, theta)

	fmt.Println("Simulation run. Now drawing images.")
	canvasWidth := 1000
	// frequency := 1000
	frequency := 100
	// scalingFactor := 1e11 // a scaling factor is needed to inflate size of stars when drawn because galaxies are very sparse
	imageList := AnimateSystem(timePoints, canvasWidth, frequency, scalingFactor)

	fmt.Println("Images drawn. Now generating GIF.")
	gifhelper.ImagesToGIF(imageList, "galaxy")
	fmt.Println("GIF drawn.")
}
