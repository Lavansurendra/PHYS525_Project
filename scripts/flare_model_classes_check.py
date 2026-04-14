import flare.model as model

print("\n--- 2. FLARE Model Classes ---")
# This prints out every specific loading tool FLARE has built-in
for item in dir(model):
    if not item.startswith("_"): # Hide the background python junk
        print(item)