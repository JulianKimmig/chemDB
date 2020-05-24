import io
import pandas as pd



class NpCharacterTool():
    long_name = None

    def __str__(self):
        return self.__class__.__name__

    def read_batch_data(self, data,name,short_name):
        raise NotImplementedError()


class MZetaNano1(NpCharacterTool):
    long_name = "Malven Zetasizer Nano"

    def read_byte_object_to_df(self, data,delimter=",",decimal="."):
        df = pd.read_csv(io.StringIO(data.decode('ansi')),delimiter=delimter,decimal=decimal)
        #rename
        import dateparser
        df.rename(columns={"Measurement Date and Time": "run_date",
                           "Sample Name": "sample_name",
                           #                          "Temperature (°C)": "temperature",
                           "Z-Average (d.nm)":"z_average",
                           "PdI":"pdi",
                           "Mean Count Rate (kcps)":"mean_count_rate",
                           "Volume Mean (d.nm)":"mean_diameter_by_volume",
                           "Intensity Mean (d.nm)":"mean_diameter_by_number",
                           "Number Mean (d.nm)":"mean_diameter_by_intensity",
                           },inplace=True)
        if "z_average" not in df.columns:
            if delimter == ",":
                return self.read_byte_object_to_df(data,delimter="\t",decimal=decimal)
            if delimter == "\t":
                return None
        if decimal == "." and "," in df["z_average"].iloc[0]:
            return self.read_byte_object_to_df(data,delimter=delimter, decimal = ",")

        df["run_date"] = df["run_date"].apply( lambda ts: dateparser.parse(ts))
        df["z_average"] = df["z_average"].astype(float).astype(int)
        #        df["temperature"] = df["temperature"].astype(float).astype(int)
        df["pdi"] = df["pdi"].astype(float)
        df["mean_count_rate"] = df["mean_count_rate"].astype(float).astype(int)
        df["mean_diameter_by_volume"] = df["mean_diameter_by_volume"].astype(float).astype(int)
        df["mean_diameter_by_number"] = df["mean_diameter_by_number"].astype(float).astype(int)
        df["mean_diameter_by_intensity"] = df["mean_diameter_by_intensity"].astype(float).astype(int)


        return df

    def read_batch_data(self, data,name,short_name):
        df = self.read_byte_object_to_df(data)

        from experiments_nanoparticle.models import NanoparticleCharacterization
        from experiments_nanoparticle.models import Nanoparticle

        batch_experiment=NanoparticleCharacterization(tool=str(self),name=name,short_name=short_name)
        sub_experiments=[]
        df = df.reset_index(drop=True)

        np=[]

        for row,data in df.iterrows():
            character = NanoparticleCharacterization(
                name="{} in {}".format(data['sample_name'], batch_experiment.short_name),
                short_name=data['sample_name'],
                tool=batch_experiment.tool,
                run_date=data['run_date'],
                batch_experiment=batch_experiment,
                batch_experiment_index=row,
            )
            sub_experiments.append(character)
            np.append(Nanoparticle(
                name=data['sample_name'],
                code=data['sample_name'].replace("  "," ").replace(" ","_"),
                z_average = data['z_average'],
                mean_diameter_by_volume = data['mean_diameter_by_volume'],
                mean_diameter_by_number = data['mean_diameter_by_number'],
                mean_diameter_by_intensity = data['mean_diameter_by_intensity'],
                pdi = data['pdi'],
#                characterizations=[character]
            ))

        return df

AVAILABLE_CHARACTERIZATIONS = {
    c.__class__.__name__: c for c in [MZetaNano1()]
}


class NanoparticleCharacterizationTool:
    choices = (
        (k, v.long_name)
        for k, v in AVAILABLE_CHARACTERIZATIONS.items()
    )
    max_length = 1


for k, v in AVAILABLE_CHARACTERIZATIONS.items():
    NanoparticleCharacterizationTool.max_length = max(NanoparticleCharacterizationTool.max_length, len(k))
    setattr(NanoparticleCharacterizationTool, k, v)

NanoparticleCharacterizationTool.max_length = 11
# print(NanoparticleCharacterizationTool.MALVEN_ZETASIZER_NANO.)
