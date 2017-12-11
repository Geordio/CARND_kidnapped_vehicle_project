/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  // TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
  //   x, y, theta and their uncertainties from GPS) and all weights to 1.
  // Add random Gaussian noise to each particle.
  // NOTE: Consult particle_filter.h for more information about this method (and others in this file).

  cout << "init" << endl;
  std::random_device rd;
  default_random_engine gen(rd());


  // tried a few values. 50 results in a higher errir. no notcable difference between 80 & 100
  // 200 performs minutely better, but has limited benefits in terms of erro (circa 0.04 for x and y
  num_particles = 200;

  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);

  // iterate though the particles
  for (int i = 0; i < num_particles; i++) {
    Particle p;
    p.id = i;
    p.x = dist_x(gen);
    p.y = dist_y(gen);
    p.theta = dist_theta(gen);
    p.weight = 1;
    particles.push_back(p);
  }

  // set the initialised flag
  is_initialized = true;
}



void ParticleFilter::prediction(double delta_t, double std_pos[],
    double velocity, double yaw_rate) {
  // TODO: Add measurements to each particle and add random Gaussian noise.
  // NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
  //  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
  //  http://www.cplusplus.com/reference/random/default_random_engine/

  cout << "prediction" << endl;
  default_random_engine gen;

  normal_distribution<double> dist_x(0, std_pos[0]);
  normal_distribution<double> dist_y(0, std_pos[1]);
  normal_distribution<double> dist_theta(0, std_pos[2]);

  // iterate though the particles
  for (int i = 0; i < num_particles; i++) {

    if (fabs(yaw_rate) < 0.00001) {
      particles[i].x += velocity * delta_t * cos(particles[i].theta);
      particles[i].y += velocity * delta_t * sin(particles[i].theta);
    } else {
      particles[i].x += velocity / yaw_rate  * (sin(particles[i].theta + yaw_rate * delta_t) - sin(particles[i].theta));
      particles[i].y += velocity / yaw_rate  * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate * delta_t));
      particles[i].theta += yaw_rate * delta_t;
    }

    // add gaussian noise
    particles[i].x += dist_x(gen);
    particles[i].y += dist_y(gen);
    particles[i].theta += dist_theta(gen);

  }

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted,
    std::vector<LandmarkObs>& observations) {
  // TODO: Find the predicted measurement that is closest to each observed measurement and assign the
  //   observed measurement to this particular landmark.
  // NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
  //   implement this method and use it as a helper during the updateWeights phase.
  //  cout << "dataAssociation" << endl;

  for (int i = 0; i < observations.size(); i++) {

    // current observations
    LandmarkObs obs = observations[i];

    // initialise the minimum distance to infinity
    double min_dist = INFINITY;
    int min_id = -1;

    for (int j = 0; j < predicted.size(); j++) {
      LandmarkObs pred = predicted[j];

      double delta_x = pred.x - obs.x;
      double delta_y = pred.y - obs.y;

      double dist_err = sqrt(delta_x * delta_x + delta_y * delta_y);

      // if this prediction is nearest to the observed landmark then set the id
      if (dist_err < min_dist) {
        min_dist = dist_err;
        min_id = pred.id;
      }
    }

    observations[i].id = min_id;

  }

}

// update weights
// Inputs:
//      sensor_range: the linit of the sensor_range
//      std_landmark: the landmark measurement uncertainty
//      observations: a vector of landmark measurements
//      map_landmarks: the map landmarks. i.e. the landmarks from the map
void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
    const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
  // TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
  //   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
  // NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
  //   according to the MAP'S coordinate system. You will need to transform between the two systems.
  //   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
  //   The following is a good resource for the theory:
  //   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
  //   and the following is a good resource for the actual equation to implement (look at equation
  //   3.33
  //   http://planning.cs.uiuc.edu/node99.html

  //  Predict the measurements to each landmark in sensor range
  // associate the data measurement to map landmarks

  cout << "sensor_range: " << sensor_range<< endl;
  cout << "observations size: " << observations.size()<< endl;
  cout << "map_landmarks size: " << map_landmarks.landmark_list.size()<< endl;

  //  cin.get();

  for (int i = 0; i < num_particles; i++) {
    double px = particles[i].x;
    double py = particles[i].y;
    double ptheta = particles[i].theta;


    vector<LandmarkObs> trans_observations;
    vector<LandmarkObs> landmk_in_range;


    // iterate over the observations and translate them to map coordinates
    for (int k = 0; k < observations.size(); k++) {

      double tx = cos(ptheta) * observations[k].x - sin(ptheta) * observations[k].y + px;
      double ty = sin(ptheta) * observations[k].x  + cos(ptheta) * observations[k].y + py;
      LandmarkObs trans_observation = { k, tx, ty };
      trans_observations.push_back(trans_observation);
    }


    // iterate over the landmarks
    // select the landmarks that are in sensor range nd stor in vector
    for (int j = 0; j < map_landmarks.landmark_list.size(); j++) {
      //      cout << "j: " << j << endl;
      // TODO get the positions from the map
      // note the positions here are floats
      float lx = map_landmarks.landmark_list[j].x_f;
      float ly = map_landmarks.landmark_list[j].y_f;
      int lid = map_landmarks.landmark_list[j].id_i;

      // calculate the difference
      float dx = lx - px;
      float dy = ly - py;
      double rmse = sqrt(dx * dx + dy * dy);

      // check if its in the sensor range, if it is then use it. If not ignore
      if ( rmse <= sensor_range) {
        LandmarkObs lm_in_range = { lid, lx, ly };
        landmk_in_range.push_back(lm_in_range);
      }
    }
    // now got transformed observations (map coords), and list of landmarks in range

    // call the data association method
    dataAssociation(landmk_in_range, trans_observations);

    //          initialise particle weight
    particles[i].weight = 1.0;

    // update the partice weights
    // iterate through the translated observations
    // then calculate the weight, based on the python code from the lessons
    // python code for lesson
    //# calculate normalization term
    //gauss_norm= (1/(2 * np.pi * sig_x * sig_y))
    //
    //# calculate exponent
    //exponent= ((x_obs - mu_x)**2)/(2 * sig_x**2) + ((y_obs - mu_y)**2)/(2 * sig_y**2)
    //
    //# calculate weight using normalization terms and exponent
    //weight= gauss_norm * math.exp(-exponent)
    //
    for (int m = 0; m < trans_observations.size(); m++) {
      double observation_id = trans_observations[m].id;
      double observation_x = trans_observations[m].x;
      double observation_y = trans_observations[m].y;

      double landmk_x;
      double landmk_y;

      double sig_x = std_landmark[0];
      double sig_y = std_landmark[1];

      // iterate through the landmarks, until we find mathcing id. probably a better way to do this
      for (unsigned int k = 0; k < landmk_in_range.size(); k++) {
        if (landmk_in_range[k].id == observation_id) {
          landmk_x = landmk_in_range[k].x;
          landmk_y = landmk_in_range[k].y;
          break;
        }
      }

      // calculate the delts from landmark to observation
      double dx = observation_x - landmk_x;
      double dy = observation_y - landmk_y;


      // calculate the weights
      double gauss_norm = (1 / (2 * M_PI * sig_x * sig_y));
      double exponent = (dx * dx) / (2 * sig_x * sig_x)  + (dy * dy) / (2 * sig_y * sig_y);

      particles[i].weight *= gauss_norm * exp(-exponent);

    }
  }
}



void ParticleFilter::resample() {
  // TODO: Resample particles with replacement with probability proportional to their weight.
  // NOTE: You may find std::discrete_distribution helpful here.
  //   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution



  // originally based on the python implementation from sebastian lesson
  //
  //  p3 = []
  //  index = int(random.random() * N)
  //  beta = 0.0
  //  mw = max(w)
  //  for i in range(N):
  //    beta += random.random() * 2.0 * mw
  //    while beta > w[index]:
  //      beta -= w[index]
  //      index = (index + 1) % N
  //    p3.append(p[index])
  //  p = p3

  std::random_device rd;
  default_random_engine gen(rd());

  vector<Particle> resampled_particles;
  vector<double> weights;
  double beta = 0;

  // get the weights from the particles
  for (int i = 0; i < num_particles; i++) {
    weights.push_back(particles[i].weight);
  }

  // used method from https://stackoverflow.com/questions/9874802/how-can-i-get-the-max-or-min-value-in-a-vector
  double max_weight = *max_element(weights.begin(), weights.end());


  // generate 2 distributions.
  // 1st of integers for the initial starting index.
  // 2nd to be used for the amount to move round the wheel
  // generate random starting index for resampling wheel using the discrete_distribution
  discrete_distribution<> dis_dist(0, num_particles);
  int index = dis_dist(gen);
  // uniform random distribution between 0 an max weight
  uniform_real_distribution<double> uni_real_dist(0.0, 1.0);

  // iterate through the particles.
  // attempt to implement the python code form the lesson
  for (int i = 0; i < num_particles; i++) {
    beta += uni_real_dist(gen) * 2.0 * max_weight;
    while (beta > weights[index]){
      beta -= weights[index];
      index = (index + 1) % num_particles;
    }
    resampled_particles.push_back(particles[index]);
  }
  particles = resampled_particles;
}


Particle ParticleFilter::SetAssociations(Particle& particle,
    const std::vector<int>& associations, const std::vector<double>& sense_x,
    const std::vector<double>& sense_y) {
  //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates

  particle.associations = associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  stringstream ss;
  copy(v.begin(), v.end(), ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length() - 1); // get rid of the trailing space
  return s;
}
string ParticleFilter::getSenseX(Particle best) {
  vector<double> v = best.sense_x;
  stringstream ss;
  copy(v.begin(), v.end(), ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length() - 1); // get rid of the trailing space
  return s;
}
string ParticleFilter::getSenseY(Particle best) {
  vector<double> v = best.sense_y;
  stringstream ss;
  copy(v.begin(), v.end(), ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length() - 1); // get rid of the trailing space
  return s;
}
