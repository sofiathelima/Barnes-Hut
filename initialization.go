package main

import (
	"math"
	"math/rand"
)

const G = 6.67408e-11 // gravitational constant -- don't change this!

const solarMass = 1.989e30 // mass of sun -- don't change this!

const blackHoleMass = 8e36 // mass of black hole -- don't change!

// InitializeUniverse() sets an initial universe given a collection of galaxies and a width.
// It returns a pointer to the resulting universe.
func InitializeUniverse(galaxies []Galaxy, w float64) *Universe {
	var u Universe
	u.width = w
	u.stars = make([]*Star, 0, len(galaxies)*len(galaxies[0].g))
	for i := range galaxies {
		u.stars = append(u.stars, galaxies[i].g...)
	}
	return &u
}

// InitializeGalaxy takes number of stars in the galaxy, radius of the galaxy to be constructed,
// and center of galaxy to be constructed. Returns a spinning Galaxy object -- which is just a slice of Star pointers
func InitializeGalaxy(numOfStars int, r, x, y float64) Galaxy {
	g := Galaxy{g: make([]*Star, numOfStars)}

	for i := range g.g {
		var s Star

		// First choose distance to center of galaxy
		dist := (rand.Float64() + 1.0) / 2.0

		// multiply by factor of r
		dist *= r

		// Next choose the angle in radians to represent the rotation
		angle := rand.Float64() * 2 * math.Pi

		// convert polar coordinates to Cartesian
		s.position.x = x + dist*math.Cos(angle)
		s.position.y = y + dist*math.Sin(angle)

		// set the mass = mass of sun by default
		s.mass = solarMass

		// set the radius equal to radius of sun in m
		s.radius = 696340000

		//set the colors
		s.red = 255
		s.green = 255
		s.blue = 255

		// now spin the galaxy

		// the following is orbital velocity equation
		//dist := Distance(pos, g[i].position)
		speed := 0.5 * math.Sqrt(G*blackHoleMass/dist) // approximation of orbital velocity equation: half of true speed to prevent instability

		s.velocity.x = speed * math.Cos(angle+math.Pi/2.0)
		s.velocity.y = speed * math.Sin(angle+math.Pi/2.0)

		//point g[i] at s
		g.g[i] = &s

	}

	//add a blackhole to the center of the galaxy

	var blackhole Star
	blackhole.mass = blackHoleMass
	blackhole.position.x = x
	blackhole.position.y = y
	blackhole.blue = 255
	blackhole.radius = 6963400000 // ten times that of a normal star (to make it visible as large)

	g.g = append(g.g, &blackhole)

	return g
}

func CreateJupiterSystem(w float64) Galaxy {
	// declaring objects
	var jupiter, io, europa, ganymede, callisto Star

	jupiter.red, jupiter.green, jupiter.blue = 223, 227, 202
	io.red, io.green, io.blue = 249, 249, 165
	europa.red, europa.green, europa.blue = 132, 83, 52
	ganymede.red, ganymede.green, ganymede.blue = 76, 0, 153
	callisto.red, callisto.green, callisto.blue = 0, 153, 76

	jupiter.mass = 1.898 * math.Pow(10, 27)
	io.mass = 8.9319 * math.Pow(10, 22)
	europa.mass = 4.7998 * math.Pow(10, 22)
	ganymede.mass = 1.4819 * math.Pow(10, 23)
	callisto.mass = 1.0759 * math.Pow(10, 23)

	jupiter.radius = 71000000
	io.radius = 1821000
	europa.radius = 1569000
	ganymede.radius = 2631000
	callisto.radius = 2410000

	jupiter.position.x, jupiter.position.y = w/2, w/2
	io.position.x, io.position.y = w/2-421600000, w/2
	europa.position.x, europa.position.y = w/2, w/2+670900000
	ganymede.position.x, ganymede.position.y = w/2+1070400000, w/2
	callisto.position.x, callisto.position.y = w/2, w/2-1882700000

	jupiter.velocity.x, jupiter.velocity.y = 0, 0
	io.velocity.x, io.velocity.y = 0, -17320
	europa.velocity.x, europa.velocity.y = -13740, 0
	ganymede.velocity.x, ganymede.velocity.y = 0, 10870
	callisto.velocity.x, callisto.velocity.y = 8200, 0

	jupiterGalaxy := Galaxy{g: []*Star{&io, &europa, &ganymede, &callisto, &jupiter}}

	return jupiterGalaxy
}
