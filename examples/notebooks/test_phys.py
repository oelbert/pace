import os
from dataclasses import dataclass, fields

import yaml

from pace import Driver, DriverConfig


config_file = "aquaplanet_c16.yaml"

with open(config_file, "r") as f:
    driver_config = DriverConfig.from_dict(yaml.safe_load(f))
driver = Driver(config=driver_config)

# for field in fields(driver.state.radiation_state):
#     try:
#         print(field.name, getattr(driver.state.radiation_state, field.name).shape)
#     except:
#         continue


try:
    driver.step_all()
finally:
    driver.cleanup()

print("done!")
