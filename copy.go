package main

func CopyUniverse(currentUniverse *Universe) *Universe {
	var newUniverse Universe

	newUniverse.width = currentUniverse.width
	newUniverse.stars = make([]*Star, len(currentUniverse.stars))

	for i := range newUniverse.stars {
		newUniverse.stars[i] = CopyStar(currentUniverse.stars[i])
	}

	return &newUniverse

}

func CopyStar(s *Star) *Star {
	var newStar Star

	// now we will make shallow copies
	newStar.position = s.position
	newStar.velocity = s.velocity
	newStar.acceleration = s.acceleration
	newStar.mass = s.mass
	newStar.radius = s.radius
	newStar.red = s.red
	newStar.green = s.green
	newStar.blue = s.blue

	return &newStar
}
