# from solcast import forecast
import numpy as np

import config

cache = {}

class Solcast:
    def __init__(
             
        self, cell_efficiency #= .17
        , area #= 4,
    ):
        self.cell_efficiency = cell_efficiency
        self.area = area

    # Function to calculate the distance between two locations (Haversine formula)
    def calculate_distance(self, lat1, lon1, lat2, lon2):
        from math import radians, sin, cos, sqrt, atan2

    # Radius of the Earth in kilometers
        earth_radius = 6371

    # Convert latitude and longitude from degrees to radians
        lat1, lon1, lat2, lon2 = map(radians, [lat1, lon1, lat2, lon2])

    # Haversine formula
        dlat = lat2 - lat1
        dlon = lon2 - lon1
        a = sin(dlat / 2) ** 2 + cos(lat1) * cos(lat2) * sin(dlon / 2) ** 2
        c = 2 * atan2(sqrt(a), sqrt(1 - a))
        distance = earth_radius * c

        return distance

    def getincident(self, latitude, longitude, globaltime):
        threshold = 0.3
        cache_key = (latitude, longitude, globaltime)
        # dict  = {}
        res = forecast.radiation_and_weather(
            hours = 15,
            api_key="9iiCFl3djay3KFYH7nWohFtkQ5lslYZE",
            latitude=latitude,
            longitude=longitude,
            output_parameters=['ghi']
            )
        a = res.to_pandas()
        # column_headings = a.columns

# Convert column headings to a list if needed
        # column_headings_list = column_headings.tolist()

# Print the column headings
        # print(column_headings)
        # print(type(a.iloc[0][0]))
        # column_headings = a['ghi'].columns
        # print(column_headings)
        # a['period_end'] = pd.to_datetime(a['ghi'])

    # # Calculate the time difference from 8 AM on the same day and convert it to minutes
        # a['minutes_elapsed'] = (a['period_end'] - a['period_end'].dt.normalize() + pd.to_timedelta('8:00:00')).dt.total_seconds() / 60

        # print(a.head)
        # print(a.head())
        intensity = a.iloc[int((globaltime//3600))]['ghi']
        print(intensity)
        # print(intensity)
        # A = [(26, 63, 1), (25, 63, 2), (24, 63, 3)]
        # A = np.array(A)
        # leftbottom = np.array((latitude, longitude, globaltime))
        # distances = np.linalg.norm(A-leftbottom, axis=1)
        # min_index = np.argmin(distances)
        # intensity = dict[A[min_index]]
        return intensity
    

# Function to retrieve data from the API or cache
    def get_data_from_api_or_cache(self,latitude, longitude, timestamp):
        threshold = 1200
        # Find the closest cache key
        closest_cache_key = None
        min_time_difference = float('inf')

        for cache_key in cache:
            cached_latitude, cached_longitude, cached_timestamp = cache_key
            location_distance = self.calculate_distance(latitude, longitude, cached_latitude, cached_longitude)
            time_difference = abs(timestamp - cached_timestamp)

            if location_distance < 0.01 and time_difference < min_time_difference:
                min_time_difference = time_difference
                closest_cache_key = cache_key

        if closest_cache_key and min_time_difference <= threshold:
            return cache[closest_cache_key]

    # If no close cache key is found or it's not close enough, make an API request
        data = self.getincident(latitude, longitude, timestamp)

    # Update the cache with the new data and timestamp
        cache[(latitude, longitude, timestamp)] = data
        return data

    def calculate_energy(self, time, globaltime, latitude, longitude):
        # Calculate power generated by solar in the path
        power = self.get_data_from_api_or_cache(latitude, longitude, globaltime)*self.cell_efficiency*self.area
        energy = power*time
        return energy


class Emulator:
    DT = config.RaceEndTime - config.RaceStartTime

    def __init__(self, cell_efficiency, area):
        self.cell_efficiency = cell_efficiency
        self.area = area

    @staticmethod
    def _calc_solar_irradiance(time):
        return 1073.099 * np.exp(-0.5 * ((time - 51908.735) / 11484.950)**2)

    def getincident(self, latitude, longitude, globaltime):
        gt = globaltime % self.DT
        intensity = self._calc_solar_irradiance(config.RaceStartTime + gt)
        return intensity
    
    # Function to retrieve data from the API or cache
    get_data_from_api_or_cache = getincident

    def calculate_energy(self, dt, globaltime, latitude, longitude):
        # Calculate power generated by solar in the path
        power = self.get_data_from_api_or_cache(latitude, longitude, globaltime) * self.cell_efficiency * self.area
        energy = power * dt / 3600
        return energy

# Solar = Solcast
Solar = Emulator